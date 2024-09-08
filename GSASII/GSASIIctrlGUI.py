# -*- coding: utf-8 -*-
#GSASIIctrlGUI - Custom GSAS-II GUI controls
'''Documentation for all the routines in module :mod:`GSASIIctrlGUI`
follows.
'''
# Documentation moved to doc/source/GSASIIGUIr.rst
#
from __future__ import division, print_function
import os
import sys
import platform
try:
    import wx
    import wx.grid as wg
    # import wx.wizard as wz
    import wx.aui
    import wx.lib.scrolledpanel as wxscroll
    import wx.html        # could postpone this for quicker startup
    import wx.lib.mixins.listctrl  as  listmix
    import wx.richtext as wxrt
    import wx.lib.filebrowsebutton as wxfilebrowse
    import matplotlib as mpl
    import matplotlib.figure as mplfig

except ImportError:
    print('ImportError for wx/mpl in GSASIIctrlGUI: ignore if docs build')
        
import time
import glob
import copy
#import ast
import random as ran
import numpy as np

import matplotlib as mpl
try:
    from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
except ImportError:
    from matplotlib.backends.backend_wx import FigureCanvas as Canvas

import GSASIIpath
import GSASIIdataGUI as G2gd
import GSASIIpwdGUI as G2pdG
import GSASIIspc as G2spc
import GSASIIobj as G2obj
import GSASIIfiles as G2fil
import GSASIIElem as G2elem
import GSASIIscriptable as G2sc
import GSASIIpwd as G2pwd
import GSASIIlattice as G2lat
import GSASIImath as G2mth
import GSASIIstrMain as G2stMn
import GSASIIIO as G2IO
import config_example
from tutorialIndex import tutorialIndex
if sys.version_info[0] >= 3:
    unicode = str
    basestring = str

# Define a short names for convenience
DULL_YELLOW = (230,230,190)
# Don't depend on wx, for scriptable
try:
    WACV = wx.ALIGN_CENTER_VERTICAL
    wxaui_NB_TOPSCROLL = wx.aui.AUI_NB_TOP | wx.aui.AUI_NB_SCROLL_BUTTONS
except:     # Don't depend on GUI
    wxaui_NB_TOPSCROLL = None

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

try:    #phoenix
    wxValidator = wx.Validator
except AttributeError:  #classic - i.e. old
    wxValidator = wx.pyValidator

#### Fixed definitions for wx Ids ################################################################################
def Define_wxId(*args):
    '''routine to create unique global wx Id symbols in this module. 
    '''
    for arg in args:
        if GSASIIpath.GetConfigValue('debug') and not arg.startswith('wxID_'):
            print (f'DBG_Problem in ID name: {arg} needs to start w/wxID_')
        if arg in globals():
            if GSASIIpath.GetConfigValue('debug'):
                print (f'DBG_Warning: {arg} already defined')
            continue
        exec(f'global {arg}; {arg} = wx.NewId()')

#### Tree Control ################################################################################
class G2TreeCtrl(wx.TreeCtrl):
    '''Create a wrapper around the standard TreeCtrl so we can "wrap"
    various events.
    
    This logs when a tree item is selected (in :meth:`onSelectionChanged`)

    This also wraps lists and dicts pulled out of the tree to track where
    they were retrieved from.
    '''
    #def SelectItem(self,event):
    #    print 'Select Item'
    #    import GSASIIobj as G2obj
    #    G2obj.HowDidIgetHere()
    #    wx.TreeCtrl.SelectItem(self,event)
        
    def __init__(self,parent=None,*args,**kwargs):
        super(self.__class__,self).__init__(parent=parent,*args,**kwargs)
        self.G2frame = parent.GetTopLevelParent()
        self.root = self.AddRoot('Loaded Data: ')
        self.SelectionChanged = None
        self.textlist = None

    def _getTreeItemsList(self,item):
        '''Get the full tree hierarchy from a reference to a tree item.
        Note that this effectively hard-codes phase and histogram names in the
        returned list. We may want to make these names relative in the future.
        '''
        textlist = [self.GetItemText(item)]
        parent = self.GetItemParent(item)
        while parent:
            if parent == self.root: break
            textlist.insert(0,self.GetItemText(parent))
            parent = self.GetItemParent(parent)
        return textlist
    
    def GetItemPyData(self,treeId):
        if 'phoenix' in wx.version():
            return wx.TreeCtrl.GetItemData(self,treeId)
        else:
            return wx.TreeCtrl.GetItemPyData(self,treeId)

    def SetItemPyData(self,treeId,data):
        if 'phoenix' in wx.version():
            return wx.TreeCtrl.SetItemData(self,treeId,data)
        else:
            return wx.TreeCtrl.SetItemPyData(self,treeId,data)

    def UpdateSelection(self):
        TId = self.GetFocusedItem()
        self.SelectItem(self.root)
        self.SelectItem(TId)

    def GetRelativeHistNum(self,histname):
        '''Returns list with a histogram type and a relative number for that
        histogram, or the original string if not a histogram
        '''
        histtype = histname.split()[0]
        if histtype != histtype.upper(): # histograms (only) have a keyword all in caps
            return histname
        item, cookie = self.GetFirstChild(self.root)
        i = 0
        while item:
            itemtext = self.GetItemText(item)
            if itemtext == histname:
                return histtype,i
            elif itemtext.split()[0] == histtype:
                i += 1
            item, cookie = self.GetNextChild(self.root, cookie)
        else:
            raise Exception("Histogram not found: "+histname)

    def ConvertRelativeHistNum(self,histtype,histnum):
        '''Converts a histogram type and relative histogram number to a
        histogram name in the current project
        '''
        item, cookie = self.GetFirstChild(self.root)
        i = 0
        while item:
            itemtext = self.GetItemText(item)
            if itemtext.split()[0] == histtype:
                if i == histnum: return itemtext
                i += 1
            item, cookie = self.GetNextChild(self.root, cookie)
        else:
            raise Exception("Histogram #'+str(histnum)+' of type "+histtype+' not found')
        
    def GetRelativePhaseNum(self,phasename):
        '''Returns a phase number if the string matches a phase name
        or else returns the original string
        '''
        item, cookie = self.GetFirstChild(self.root)
        while item:
            itemtext = self.GetItemText(item)
            if itemtext == "Phases":
                parent = item
                item, cookie = self.GetFirstChild(parent)
                i = 0
                while item:
                    itemtext = self.GetItemText(item)
                    if itemtext == phasename:
                        return i
                    item, cookie = self.GetNextChild(parent, cookie)
                    i += 1
                else:
                    return phasename # not a phase name
            item, cookie = self.GetNextChild(self.root, cookie)
        else:
            raise Exception("No phases found ")

    def ConvertRelativePhaseNum(self,phasenum):
        '''Converts relative phase number to a phase name in
        the current project
        '''
        item, cookie = self.GetFirstChild(self.root)
        while item:
            itemtext = self.GetItemText(item)
            if itemtext == "Phases":
                parent = item
                item, cookie = self.GetFirstChild(parent)
                i = 0
                while item:
                    if i == phasenum:
                        return self.GetItemText(item)
                    item, cookie = self.GetNextChild(parent, cookie)
                    i += 1
                else:
                    raise Exception("Phase "+str(phasenum)+" not found")
            item, cookie = self.GetNextChild(self.root, cookie)
        else:
            raise Exception("No phases found ")

    def GetImageLoc(self,TreeId):
        '''Get Image data from the Tree. Handles cases where the
        image name is specified, as well as where the image file name is
        a tuple containing the image file and an image number
        '''
        
        size,imagefile = self.GetItemPyData(TreeId)
        if type(imagefile) is tuple or type(imagefile) is list:
            return size,imagefile[0],imagefile[1]
        else:
            return size,imagefile,None

    def UpdateImageLoc(self,TreeId,imagefile):
        '''Saves a new imagefile name in the Tree. Handles cases where the
        image name is specified, as well as where the image file name is
        a tuple containing the image file and an image number
        '''
        
        idata = self.GetItemPyData(TreeId)
        if type(idata[1]) is tuple or type(idata[1]) is list:
            idata[1] = list(idata[1])
            idata[1][0] = [imagefile,idata[1][1]]
        else:
            idata[1]  = imagefile
        
    def SaveExposedItems(self):
        '''Traverse the top level tree items and save names of exposed (expanded) tree items.
        Done before a refinement.
        '''
        self.ExposedItems = []
        item, cookie = self.GetFirstChild(self.root)
        while item:
            name = self.GetItemText(item)
            if self.IsExpanded(item): self.ExposedItems.append(name)
            item, cookie = self.GetNextChild(self.root, cookie)
#        print 'exposed:',self.ExposedItems

    def RestoreExposedItems(self):
        '''Traverse the top level tree items and restore exposed (expanded) tree items
        back to their previous state (done after a reload of the tree after a refinement)
        '''
        item, cookie = self.GetFirstChild(self.root)
        while item:
            name = self.GetItemText(item)
            if name in self.ExposedItems: self.Expand(item)
            item, cookie = self.GetNextChild(self.root, cookie)

def ReadOnlyTextCtrl(*args,**kwargs):
    '''Create a read-only TextCtrl for display of constants
    This is probably not ideal as it mixes visual cues, but it does look nice.
    Addresses 4.2 bug where TextCtrl has no default size
    '''
    kwargs['style'] = wx.TE_READONLY
    if wx.__version__.startswith('4.2') and 'size' not in kwargs:
        kwargs['size'] = (105, 22)
    Txt = wx.TextCtrl(*args,**kwargs)
    Txt.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE))
    return Txt
            
#### TextCtrl that stores input as entered with optional validation ################################################################################
class ValidatedTxtCtrl(wx.TextCtrl):
    '''Create a TextCtrl widget that uses a validator to prevent the
    entry of inappropriate characters and changes color to highlight
    when invalid input is supplied. As valid values are typed,
    they are placed into the dict or list where the initial value
    came from. The type of the initial value must be int,
    float or str or None (see :obj:`key` and :obj:`typeHint`);
    this type (or the one in :obj:`typeHint`) is preserved.

    Float values can be entered in the TextCtrl as numbers or also
    as algebraic expressions using operators + - / \\* () and \\*\\*,
    in addition pi, sind(), cosd(), tand(), and sqrt() can be used,
    as well as appreviations s, sin, c, cos, t, tan and sq. 

    :param wx.Panel parent: name of panel or frame that will be
      the parent to the TextCtrl. Can be None.

    :param dict/list loc: the dict or list with the initial value to be
      placed in the TextCtrl. 

    :param int/str key: the dict key or the list index for the value to be
      edited by the TextCtrl. The ``loc[key]`` element must exist, but may
      have value None. If None, the type for the element is taken from
      :obj:`typeHint` and the value for the control is set initially
      blank (and thus invalid.) This is a way to specify a field without a
      default value: a user must set a valid value.
       
      If the value is not None, it must have a base
      type of int, float, str or unicode; the TextCrtl will be initialized
      from this value.
      
    :param list nDig: number of digits, places and optionally the format 
       ([nDig,nPlc,fmt]) after decimal to use for display of float. The format 
       is either 'f' (default) or 'g'. Alternately, None can be specified which 
       causes numbers to be displayed with approximately 5 significant figures
       for floats. If this is specified, then :obj:`typeHint` = float becomes the
       default. 
       (Default=None).

    :param bool notBlank: if True (default) blank values are invalid
      for str inputs.
      
    :param number xmin: minimum allowed valid value. If None (default) the
      lower limit is unbounded.
      NB: test in NumberValidator is val >= xmin not val > xmin

    :param number xmax: maximum allowed valid value. If None (default) the
      upper limit is unbounded
      NB: test in NumberValidator is val <= xmax not val < xmax
      
    :param list exclLim: if True exclude min/max value ([exclMin,exclMax]); 
      (Default=[False,False]) 

    :param function OKcontrol: specifies a function or method that will be
      called when the input is validated. The called function is supplied
      with one argument which is False if the TextCtrl contains an invalid
      value and True if the value is valid.
      Note that this function should check all values
      in the dialog when True, since other entries might be invalid.
      The default for this is None, which indicates no function should
      be called.

    :param function OnLeave: specifies a function or method that will be
      called when the focus for the control is lost.
      The called function is supplied with (at present) three keyword arguments:

      * invalid: (*bool*) True if the value for the TextCtrl is invalid
      * value:   (*int/float/str*)  the value contained in the TextCtrl
      * tc:      (*wx.TextCtrl*)  the TextCtrl object

      The number of keyword arguments may be increased in the future should needs arise,
      so it is best to code these functions with a \\*\\*kwargs argument so they will
      continue to run without errors

      The default for OnLeave is None, which indicates no function should
      be called.

    :param type typeHint: the value of typeHint should be int, float or str (or None).
      The value for this will override the initial type taken from value
      for the dict/list element ``loc[key]`` if not None and thus specifies the
      type for input to the TextCtrl.
      Defaults as None, which is ignored, unless  :obj:`nDig` is specified in which
      case the default is float.

    :param bool CIFinput: for str input, indicates that only printable
      ASCII characters may be entered into the TextCtrl. Forces output
      to be ASCII rather than Unicode. For float and int input, allows
      use of a single '?' or '.' character as valid input.

    :param dict OnLeaveArgs: a dict with keyword args that are passed to
      the :attr:`OnLeave` function. Defaults to ``{}``

    :param bool ASCIIonly: if set as True will remove unicode characters from
      strings

    :param (other): other optional keyword parameters for the
      wx.TextCtrl widget such as size or style may be specified.
    '''
    def __init__(self,parent,loc,key,nDig=None,notBlank=True,xmin=None,xmax=None,
        OKcontrol=None,OnLeave=None,typeHint=None,CIFinput=False,exclLim=[False,False],
        OnLeaveArgs={}, ASCIIonly=False,
                     min=None, max=None, # patch: remove this eventually
                     **kw):
        # save passed values needed outside __init__
        self.result = loc
        self.key = key
        self.nDig = nDig
        self.OKcontrol=OKcontrol
        self.OnLeave = OnLeave
        self.OnLeaveArgs = OnLeaveArgs
        self.CIFinput = CIFinput
        self.notBlank = notBlank
        self.ASCIIonly = ASCIIonly
        
        # patch: remove this when min & max are no longer used to call this
        if min is not None:
            xmin=min
            if GSASIIpath.GetConfigValue('debug'):
                print('Call to ValidatedTxtCtrl using min (change to xmin) here:')
                G2obj.HowDidIgetHere(True)
        if max is not None:
            xmax=max
            if GSASIIpath.GetConfigValue('debug'):
                print('Call to ValidatedTxtCtrl using max (change to xmax) here:')
                G2obj.HowDidIgetHere(True)
        # end patch
        
        # initialization
        self.invalid = False   # indicates if the control has invalid contents
        self.evaluated = False # set to True when the validator recognizes an expression
        self.timer = None      # tracks pending updates for expressions in float textctrls
        self.delay = 5000      # delay for timer update (5 sec)
        self.type = str
        
        val = loc[key]
        if 'style' in kw: # add a "Process Enter" to style
            kw['style'] |= wx.TE_PROCESS_ENTER
        else:
            kw['style'] = wx.TE_PROCESS_ENTER
        if 'size' not in kw: # wx 4.2.0 needs a size
            kw['size'] = (105,-1)
        if typeHint is not None:
            self.type = typeHint
        elif nDig is not None:
            self.type = float
        elif 'int' in str(type(val)):
            self.type = int
        elif 'float' in str(type(val)):
            self.type = float
        elif isinstance(val,str) or isinstance(val,unicode):
            self.type = str
        elif val is None:
            raise Exception("ValidatedTxtCtrl error: value of "+str(key)+
                             " element is None and typeHint not defined as int or float")
        else:
            raise Exception("ValidatedTxtCtrl error: Unknown element ("+str(key)+
                             ") type: "+str(type(val)))
        if self.type is int:        
            wx.TextCtrl.__init__(self,parent,wx.ID_ANY,
                validator=NumberValidator(int,result=loc,key=key,xmin=xmin,xmax=xmax,
                    exclLim=exclLim,OKcontrol=OKcontrol,CIFinput=CIFinput),**kw)
            if val is not None:
                self._setValue(val)
            else: # no default is invalid for a number
                self.invalid = True
                self._IndicateValidity()
        elif self.type is float:
            wx.TextCtrl.__init__(self,parent,wx.ID_ANY,
                validator=NumberValidator(float,result=loc,key=key,xmin=xmin,xmax=xmax,
                    exclLim=exclLim,OKcontrol=OKcontrol,CIFinput=CIFinput),**kw)
            if val is not None:
                self._setValue(val)
            else:
                self.invalid = True
                self._IndicateValidity()
        else:
            if self.CIFinput:
                wx.TextCtrl.__init__(
                    self,parent,wx.ID_ANY,
                    validator=ASCIIValidator(result=loc,key=key),
                    **kw)
            else:
                wx.TextCtrl.__init__(self,parent,wx.ID_ANY,**kw)
            if val is not None:
                self.SetValue(val)
            if notBlank:
                self.Bind(wx.EVT_CHAR,self._onStringKey)
                self.ShowStringValidity() # test if valid input
            else:
                self.invalid = False
                self.Bind(wx.EVT_CHAR,self._GetStringValue)
        
        # When the mouse is moved away or the widget loses focus,
        # display the last saved value, if an expression
        self.Bind(wx.EVT_LEAVE_WINDOW, self._onLeaveWindow)
        self.Bind(wx.EVT_KILL_FOCUS, self._onLoseFocus)
        self.Bind(wx.EVT_TEXT_ENTER, self._onLoseFocus)
        # patch for wx 2.9 on Mac
        i,j= wx.__version__.split('.')[0:2]
        if int(i)+int(j)/10. > 2.8 and 'wxOSX' in wx.PlatformInfo:
            self.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)

    def SetValue(self,val):
        if self.result is not None: # note that this bypasses formatting
            self.result[self.key] = val
        self._setValue(val)

    def _setValue(self,val,show=True):
        '''Check the validity of an int or float value and convert to a str.
        Possibly format it. If show is True, display the formatted value in
        the Text widget.
        '''
        self.invalid = False
        if self.type is int:
            try:
                if int(val) != val:
                    self.invalid = True
                else:
                    val = int(val)
            except:
                if self.CIFinput and (val == '?' or val == '.'):
                    pass
                else:
                    self.invalid = True
            if show and not self.invalid: wx.TextCtrl.SetValue(self,str(val))
        elif self.type is float:
            try:
                if type(val) is str: val = val.replace(',','.')
                val = float(val) # convert strings, if needed
            except:
                if self.CIFinput and (val == '?' or val == '.'):
                    pass
                else:
                    self.invalid = True
            if self.nDig and show and not self.invalid:
                wx.TextCtrl.SetValue(self,str(G2fil.FormatValue(val,self.nDig)))
                self.evaluated = False # expression has been recast as value, reset flag
            elif show and not self.invalid:
                wx.TextCtrl.SetValue(self,str(G2fil.FormatSigFigs(val)).rstrip('0'))
                self.evaluated = False # expression has been recast as value, reset flag
        else:
            if self.ASCIIonly:
                s = ''
                for c in val:
                    if ord(c) < 128:
                        s += c
                    else:
                        s += '!'
                if val != s:
                    val = s
                    show = True
            if show:
                try:
                    wx.TextCtrl.SetValue(self,str(val))
                except:
                    wx.TextCtrl.SetValue(self,val)
            self.ShowStringValidity() # test if valid input
            return
        
        self._IndicateValidity()
        if self.OKcontrol:
            self.OKcontrol(not self.invalid)

    def OnKeyDown(self,event):
        'Special callback for wx 2.9+ on Mac where backspace is not processed by validator'
        key = event.GetKeyCode()
        if key in [wx.WXK_BACK, wx.WXK_DELETE]:
            if self.Validator: wx.CallAfter(self.Validator.TestValid,self)
        if key == wx.WXK_RETURN or key == wx.WXK_NUMPAD_ENTER:
            self._onLoseFocus(None)
        if event: event.Skip()
        if self.timer:
            self.timer.Restart(self.delay)
                    
    def _onStringKey(self,event):
        if event: event.Skip()
        if self.invalid: # check for validity after processing the keystroke
            wx.CallAfter(self.ShowStringValidity,True) # was invalid
        else:
            wx.CallAfter(self.ShowStringValidity,False) # was valid

    def _IndicateValidity(self):
        'Set the control colors to show invalid input'
        if self.invalid:
            ins = self.GetInsertionPoint()
            self.SetForegroundColour("red")
            self.SetBackgroundColour("yellow")
            self.SetFocus()
            self.Refresh() # this selects text on some Linuxes
            self.SetSelection(0,0)   # unselect
            self.SetInsertionPoint(ins) # put insertion point back 
        else: # valid input
            self.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            self.SetForegroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT))
            self.Refresh()
            self.SetFocus() # seems needed, at least on MacOS to get color change

    def _GetNumValue(self):
        'Get and where needed convert string from GetValue into int or float'
        try:
            val = self.GetValue()
            if self.type is int:
                val = int(val)
            elif self.type is float:
                val = float(val)
        except:
            if self.CIFinput and (val == '?' or val == '.'):
                pass
            else:
                self.invalid = True
        return val
    
    def ShowStringValidity(self,previousInvalid=True):
        '''Check if input is valid. Anytime the input is
        invalid, call self.OKcontrol (if defined) because it is fast.
        If valid, check for any other invalid entries only when
        changing from invalid to valid, since that is slower.
        
        :param bool previousInvalid: True if the TextCtrl contents were
          invalid prior to the current change.
          
        '''
        val = self.GetValue().strip()
        if self.notBlank:
            self.invalid = not val
        else:
            self.invalid = False
        self._IndicateValidity()
        if self.invalid:
            if self.OKcontrol:
                self.OKcontrol(False)
        elif self.OKcontrol and previousInvalid:
            self.OKcontrol(True)
        self._SaveStringValue()         # always store the result

    def _GetStringValue(self,event):
        '''Get string input and store.
        '''
        if event: event.Skip() # process keystroke
        wx.CallAfter(self._SaveStringValue)
        
    def _SaveStringValue(self):
        val = self.GetValue().strip()
        # always store the result
        if self.CIFinput and '2' in platform.python_version_tuple()[0]: # Py2/CIF make results ASCII
            self.result[self.key] = val.encode('ascii','replace') 
        else:
            self.result[self.key] = val

    def _onLeaveWindow(self,event):
        '''If the mouse leaves the text box, save the result, if valid,
        but (unlike _onLoseFocus) there is a two second delay before 
        the textbox contents are updated with the value from the formula.
        '''
        def delayedUpdate():
            self.timer = None
            try:
                self._setValue(self.result[self.key])
            except:
                pass
        if self.type is not str:
            if not self.IsModified(): return  #ignore mouse crusing
        elif self.result[self.key] == self.GetValue(): # .IsModified() seems unreliable for str 
           return
        if self.evaluated and not self.invalid: # deal with computed expressions
            if self.timer:
                self.timer.Restart(self.delay)
            else:
                self.timer = wx.CallLater(self.delay,delayedUpdate)    # this includes a try in case the widget is deleted
        if self.invalid: # don't update an invalid expression
            if event: event.Skip()
            return
        self._setValue(self.result[self.key],show=False) # save value quietly
        if self.OnLeave:
            self.event = event
            self.OnLeave(invalid=self.invalid,value=self.result[self.key],
                tc=self,**self.OnLeaveArgs)
        if event: event.Skip()
            
    def _onLoseFocus(self,event):
        '''Enter has been pressed or focus transferred to another control,
        Evaluate and update the current control contents
        '''
        if event: event.Skip()
        if self.type is not str:
            if not self.IsModified(): return  #ignore mouse crusing
        elif self.result[self.key] == self.GetValue(): # .IsModified() seems unreliable for str 
           return
        if self.evaluated: # deal with computed expressions
            if self.invalid: # don't substitute for an invalid expression
                return 
            self._setValue(self.result[self.key])
        elif self.result is not None: # show formatted result, as Bob wants
            self.result[self.key] = self._GetNumValue()
            if not self.invalid: # don't update an invalid expression
                self._setValue(self.result[self.key])
                
        if self.OnLeave:
            self.event = event
            try:
                self.OnLeave(invalid=self.invalid,value=self.result[self.key],
                    tc=self,**self.OnLeaveArgs)
            except:
                pass
################################################################################
class NumberValidator(wxValidator):
    '''A validator to be used with a TextCtrl to prevent
    entering characters other than digits, signs, and for float
    input, a period and exponents.
    
    The value is checked for validity after every keystroke
      If an invalid number is entered, the box is highlighted.
      If the number is valid, it is saved in result[key]

    :param type typ: the base data type. Must be int or float.

    :param bool positiveonly: If True, negative integers are not allowed
      (default False). This prevents the + or - keys from being pressed.
      Used with typ=int; ignored for typ=float.

    :param number xmin: Minimum allowed value. If None (default) the
      lower limit is unbounded

    :param number xmax: Maximum allowed value. If None (default) the
      upper limit is unbounded
      
    :param list exclLim: if True exclude xmin/xmax value ([exclMin,exclMax]); 
     (Default=[False,False]) 

    :param dict/list result: List or dict where value should be placed when valid

    :param any key: key to use for result (int for list)

    :param function OKcontrol: function or class method to control
      an OK button for a window. 
      Ignored if None (default)

    :param bool CIFinput: allows use of a single '?' or '.' character
      as valid input.
      
    '''
    def __init__(self, typ, positiveonly=False, xmin=None, xmax=None,exclLim=[False,False],
        result=None, key=None, OKcontrol=None, CIFinput=False):
        'Create the validator'
        wxValidator.__init__(self)
        # save passed parameters
        self.typ = typ
        self.positiveonly = positiveonly
        self.xmin = xmin
        self.xmax = xmax
        self.exclLim = exclLim
        self.result = result
        self.key = key
        self.OKcontrol = OKcontrol
        self.CIFinput = CIFinput
        # set allowed keys by data type
        self.Bind(wx.EVT_CHAR, self.OnChar)
        if self.typ == int and self.positiveonly:
            self.validchars = '0123456789'
        elif self.typ == int:
            self.validchars = '0123456789+-'
        elif self.typ == float:
            # allow for above and sind, cosd, sqrt, tand, pi, and abbreviations
            # also addition, subtraction, division, multiplication, exponentiation
            self.validchars = '0123456789.-+eE/cosindcqrtap()*,'
        else:
            self.validchars = None
            return
        if self.CIFinput:
            self.validchars += '?.'
    def Clone(self):
        'Create a copy of the validator, a strange, but required component'
        return NumberValidator(typ=self.typ, 
                               positiveonly=self.positiveonly,
                               xmin=self.xmin, xmax=self.xmax,
                               result=self.result, key=self.key,
                               OKcontrol=self.OKcontrol,
                               CIFinput=self.CIFinput)
    def TransferToWindow(self):
        'Needed by validator, strange, but required component'
        return True # Prevent wxDialog from complaining.
    def TransferFromWindow(self):
        'Needed by validator, strange, but required component'
        return True # Prevent wxDialog from complaining.
    def TestValid(self,tc):
        '''Check if the value is valid by casting the input string
        into the current type.

        Set the invalid variable in the TextCtrl object accordingly.

        If the value is valid, save it in the dict/list where
        the initial value was stored, if appropriate. 

        :param wx.TextCtrl tc: A reference to the TextCtrl that the validator
          is associated with.
        '''
        tc.invalid = False # assume valid
        if self.CIFinput:
            val = tc.GetValue().strip()
            if val == '?' or val == '.':
                self.result[self.key] = val
                return
        try:
            val = self.typ(tc.GetValue())
        except (ValueError, SyntaxError):
            if self.typ is float: # for float values, see if an expression can be evaluated
                val = G2fil.FormulaEval(tc.GetValue().replace(',','.'))
                if val is None:
                    tc.invalid = True
                    return
                else:
                    tc.evaluated = True
            else: 
                tc.invalid = True
                return
        if self.xmax != None:
            if val >= self.xmax and self.exclLim[1]:
                tc.invalid = True
            elif val > self.xmax:
                tc.invalid = True
        if self.xmin != None:
            if val <= self.xmin and self.exclLim[0]:
                tc.invalid = True
            elif val < self.xmin:
                tc.invalid = True  # invalid
        if self.key is not None and self.result is not None and not tc.invalid:
            self.result[self.key] = val

    def ShowValidity(self,tc):
        '''Set the control colors to show invalid input

        :param wx.TextCtrl tc: A reference to the TextCtrl that the validator
          is associated with.

        '''
        if tc.invalid:
            ins = tc.GetInsertionPoint()
            tc.SetForegroundColour("red")
            tc.SetBackgroundColour("yellow")
            tc.SetFocus()
            tc.Refresh() # this selects text on some Linuxes
            tc.SetSelection(0,0)   # unselect
            tc.SetInsertionPoint(ins) # put insertion point back 
            return False
        else: # valid input
            tc.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            tc.SetForegroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT))
            tc.Refresh()
            return True

    def CheckInput(self,previousInvalid):
        '''called to test every change to the TextCtrl for validity and
        to change the appearance of the TextCtrl

        Anytime the input is invalid, call self.OKcontrol
        (if defined) because it is fast. 
        If valid, check for any other invalid entries only when
        changing from invalid to valid, since that is slower.

        :param bool previousInvalid: True if the TextCtrl contents were
          invalid prior to the current change.
        '''
        tc = self.GetWindow()
        self.TestValid(tc)
        self.ShowValidity(tc)
        # if invalid
        if tc.invalid and self.OKcontrol:
            self.OKcontrol(False)
        if not tc.invalid and self.OKcontrol and previousInvalid:
            self.OKcontrol(True)

    def OnChar(self, event):
        '''Called each type a key is pressed
        ignores keys that are not allowed for int and float types
        '''
        key = event.GetKeyCode()
        tc = self.GetWindow()
        if key == wx.WXK_RETURN or key == wx.WXK_NUMPAD_ENTER:
            if tc.invalid:
                self.CheckInput(True) 
            else:
                self.CheckInput(False) 
            if event: event.Skip()
            return
        if key < wx.WXK_SPACE or key == wx.WXK_DELETE or key > 255: # control characters get processed
            if event: event.Skip()
            if tc.invalid:
                wx.CallAfter(self.CheckInput,True) 
            else:
                wx.CallAfter(self.CheckInput,False) 
            return
        elif chr(key) in self.validchars: # valid char pressed?
            if event: event.Skip()
            if tc.invalid:
                wx.CallAfter(self.CheckInput,True) 
            else:
                wx.CallAfter(self.CheckInput,False) 
            return
        return  # Returning without calling event.Skip, which eats the keystroke

################################################################################
class ASCIIValidator(wxValidator):
    '''A validator to be used with a TextCtrl to prevent
    entering characters other than ASCII characters.
    
    The value is checked for validity after every keystroke
      If an invalid number is entered, the box is highlighted.
      If the number is valid, it is saved in result[key]

    :param dict/list result: List or dict where value should be placed when valid

    :param any key: key to use for result (int for list)

    '''
    def __init__(self, result=None, key=None):
        'Create the validator'
        import string
        wxValidator.__init__(self)
        # save passed parameters
        self.result = result
        self.key = key
        self.validchars = string.ascii_letters + string.digits + string.punctuation + string.whitespace
        self.Bind(wx.EVT_CHAR, self.OnChar)
    def Clone(self):
        'Create a copy of the validator, a strange, but required component'
        return ASCIIValidator(result=self.result, key=self.key)
        tc = self.GetWindow()
        tc.invalid = False # make sure the validity flag is defined in parent
    def TransferToWindow(self):
        'Needed by validator, strange, but required component'
        return True # Prevent wxDialog from complaining.
    def TransferFromWindow(self):
        'Needed by validator, strange, but required component'
        return True # Prevent wxDialog from complaining.
    def TestValid(self,tc):
        '''Check if the value is valid by casting the input string
        into ASCII. 

        Save it in the dict/list where the initial value was stored

        :param wx.TextCtrl tc: A reference to the TextCtrl that the validator
          is associated with.
        '''
        if '2' in platform.python_version_tuple()[0]:
            self.result[self.key] = tc.GetValue().encode('ascii','replace')
        else:
            self.result[self.key] = tc.GetValue()

    def OnChar(self, event):
        '''Called each type a key is pressed
        ignores keys that are not allowed for int and float types
        '''
        key = event.GetKeyCode()
        tc = self.GetWindow()
        if key == wx.WXK_RETURN or key == wx.WXK_NUMPAD_ENTER:
            self.TestValid(tc)
            if event: event.Skip()
            return
        if key < wx.WXK_SPACE or key == wx.WXK_DELETE or key > 255: # control characters get processed
            if event: event.Skip()
            self.TestValid(tc)
            return
        elif chr(key) in self.validchars: # valid char pressed?
            if event: event.Skip()
            self.TestValid(tc)
            return
        return  # Returning without calling event.Skip, which eats the keystroke

class G2Slider(wx.Slider):
    '''Wrapper around wx.Slider widget that implements scaling
    Also casts floats as integers to avoid py3.10+ errors
    '''
    global ci
    ci = lambda x: int(x + 0.5)  # closest integer       
    def __init__(self, parent, id=wx.ID_ANY, value=0, minValue=0, maxValue=100, *arg, **kwarg):
        wx.Slider.__init__(self, parent, id, ci(value), ci(minValue), ci(maxValue), *arg, **kwarg)
        self.iscale = 1
        
    def SetScaling(self,iscale):
        self.iscale = iscale
        
    def SetScaledRange(self,xmin,xmax):
        self.SetRange(ci(xmin*self.iscale),ci(xmax*self.iscale))

    def SetScaledValue(self,value):
        wx.Slider.SetValue(self, ci(self.iscale*value))

    def GetScaledValue(self):
        return self.GetValue()/float(self.iscale)
    
    def SetValue(self,value):
        wx.Slider.SetValue(self, ci(value))
        
    def SetMax(self,xmax):
        wx.Slider.SetMax(self,ci(xmax*self.iscale))
        
    def SetMin(self,xmin):
        wx.Slider.SetMin(self,ci(xmin*self.iscale))

def G2SliderWidget(parent,loc,key,label,xmin,xmax,iscale,
    onChange=None,onChangeArgs=[],sizer=None,nDig=None,size=(50,20)):
    '''A customized combination of a wx.Slider and a validated 
    wx.TextCtrl (see :class:`ValidatedTxtCtrl`) that allows either
    a slider or text entry to set a value within a range.

    :param wx.Panel parent: name of panel or frame that will be
      the parent to the TextCtrl. Can be None.

    :param dict/list loc: the dict or list with the initial value to be
      placed in the TextCtrl.

    :param int/str key: the dict key or the list index for the value to be
      edited by the TextCtrl. The ``loc[key]`` element must exist and should 
      have a float value. It will be forced to an initial value 
      between xmin and xmax.
      
    :param str label: A label to be placed to the left of the slider. 

    :param float xmin: the minimum allowed valid value. 

    :param float xmax: the maximum allowed valid value.

    :param float iscale: number to scale values to integers, which is what the 
       Scale widget uses. If the xmin=1 and xmax=4 and iscale=1 then values
       only the values 1,2,3 and 4 can be set with the slider. However, 
       if iscale=2 then the values 1, 1.5, 2, 2.5, 3, 3.5 and 4 are all allowed.

    :param callable onChange: function to call when value is changed.
       Default is None where nothing will be called. 
       
    :param list onChangeArgs: arguments to be passed to onChange function 
       when called.
    :returns: returns a wx.BoxSizer containing the widgets
    '''

    def onScale(event):
        loc[key] = vScale.GetScaledValue()
        wx.TextCtrl.SetValue(vEntry,str(loc[key])) # will not trigger onValSet
        if onChange: onChange(*onChangeArgs)
    def onValSet(*args,**kwargs):
        vScale.SetScaledValue(loc[key])
        if onChange: onChange(*onChangeArgs)
    loc[key] = min(xmax,max(xmin,loc[key]))
    if sizer is None:
        hSizer = wx.BoxSizer(wx.HORIZONTAL)
    else:
        hSizer = sizer
    hSizer.Add(wx.StaticText(parent,wx.ID_ANY,label),0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    vScale = G2Slider(parent,style=wx.SL_HORIZONTAL,size=(200,25))
    vScale.SetScaling(iscale)
    vScale.SetScaledRange(xmin,xmax)
    vScale.SetScaledValue(loc[key])
    vScale.Bind(wx.EVT_SLIDER, onScale)
    if nDig is None:
        nDig = (10,int(0.9+np.log10(iscale)))
    vEntry = ValidatedTxtCtrl(parent,loc,key,nDig=nDig,OnLeave=onValSet,
                xmin=xmin,xmax=xmax,typeHint=float,size=size)
    if sizer is None:
        hSizer.Add(vEntry,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL,5)
        hSizer.Add(vScale,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL,0)
        return hSizer
    else:
        hSizer.Add(vEntry,0,wx.RIGHT|wx.ALIGN_CENTER_VERTICAL,5)
        hSizer.Add(vScale,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
        return vEntry,vScale
    
def G2SpinWidget(parent,loc,key,label,xmin=None,xmax=None,
        onChange=None,onChangeArgs=[],hsize=35):
    '''A customized combination of a wx.SpinButton and a validated 
    wx.TextCtrl (see :class:`ValidatedTxtCtrl`) that allows either
    a the spin button or text entry to set a value within a range.

    :param wx.Panel parent: name of panel or frame that will be
      the parent to the TextCtrl. Can be None.

    :param dict/list loc: the dict or list with the initial value to be
      placed in the TextCtrl.

    :param int/str key: the dict key or the list index for the value to be
      edited by the TextCtrl. The ``loc[key]`` element must exist and should 
      have a float or int value. It will be forced to an integer initial value 
      between xmin and xmax.
      
    :param str label: A label to be placed to the left of the entry widget. 

    :param int xmin: the minimum allowed valid value. If None it is ignored.

    :param int xmax: the maximum allowed valid value. If None it is ignored.

    :param callable onChange: function to call when value is changed.
       Default is None where nothing will be called. 
       
    :param list onChangeArgs: arguments to be passed to onChange function 
       when called.

    :param int hsize: length of TextCtrl in pixels. Defaults to 35.

    :returns: returns a wx.BoxSizer containing the widgets
    '''

    def _onSpin(event):
        Obj = event.GetEventObject()
        loc[key] += Obj.GetValue()  # +1 or -1
        if xmin is not None and loc[key] < xmin:
            loc[key] = xmin
        if xmax is not None and loc[key] > xmax:
            loc[key] = xmax
        wx.TextCtrl.SetValue(vEntry,str(loc[key])) # will not trigger onValSet
        Obj.SetValue(0)
        if onChange: onChange(*onChangeArgs)
    def _onValSet(*args,**kwargs):
        if onChange: onChange(*onChangeArgs)
    if xmin is not None:
        loc[key] = max(xmin,loc[key])
    if xmax is not None:
        loc[key] = min(xmax,loc[key])
    hSizer = wx.BoxSizer(wx.HORIZONTAL)
    if label: 
        hSizer.Add(wx.StaticText(parent,wx.ID_ANY,label),0,
                       wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    spin = wx.SpinButton(parent,style=wx.SP_VERTICAL,size=wx.Size(20,20))
    spin.SetRange(-1,1)
    spin.Bind(wx.EVT_SPIN, _onSpin)
    loc[key] = int(loc[key]+0.5)
    vEntry = ValidatedTxtCtrl(parent,loc,key,OnLeave=_onValSet,
                xmin=xmin,xmax=xmax,typeHint=int,size=(hsize,-1))
    hSizer.Add(vEntry,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL,5)
    hSizer.Add(spin,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    return hSizer

################################################################################
def HorizontalLine(sizer,parent):
    '''Draws a horizontal line as wide as the window.
    '''
    if sys.platform == "darwin": 
        #sizer.Add((-1,2))
        line = wx.Panel(parent, size=(-1, 2))
        #line.SetBackgroundColour('red')
        line.SetBackgroundColour((128,128,128))
        sizer.Add(line, 0, wx.EXPAND|wx.ALL, 0)
        sizer.Add((-1,5))
    else:
        line = wx.StaticLine(parent, size=(-1,3), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.EXPAND|wx.ALL, 5)

################################################################################
class G2Button(wx.Button):
    '''A version of wx.Button. Bindings are saved
    in the object, and are looked up rather than directly set with a bind.
    :param wx.Panel parent: parent widget
    :param int id: Id for button
    :param str label: label for button
    :param function handler: a routine to call when the button is pressed
    '''
    def __init__(self,parent,id=wx.ID_ANY,label='',locationcode='',
                 handler=None,*args,**kwargs):
        super(self.__class__,self).__init__(parent,id,label,*args,**kwargs)
        self.label = label
        self.handler = handler
        self.locationcode = locationcode
        self.Bind(wx.EVT_BUTTON,self.onPress)
    def onPress(self,event):
        'create log event and call handler'
        self.handler(event)
        
################################################################################
class EnumSelector(wx.ComboBox):
    '''A customized :class:`wxpython.ComboBox` that selects items from a list
    of choices, but sets a dict (list) entry to the corresponding
    entry from the input list of values.

    :param wx.Panel parent: the parent to the :class:`~wxpython.ComboBox` (usually a
      frame or panel)
    :param dct: a dict or list to contain the value set
      for the :class:`~wxpython.ComboBox`.
    :param item: the dict key (or list index) where ``dct[item]`` will 
      be set to the value selected in the :class:`~wxpython.ComboBox`. Also, dct[item]
      contains the starting value shown in the widget. If the value
      does not match an entry in :data:`values`, the first value
      in :data:`choices` is used as the default, but ``dct[item]`` is
      not changed.    
    :param list choices: a list of choices to be displayed to the
      user such as
      ::
      
      ["default","option 1","option 2",]

      Note that these options will correspond to the entries in 
      :data:`values` (if specified) item by item. 
    :param list values: a list of values that correspond to
      the options in :data:`choices`, such as
      ::
      
      [0,1,2]
      
      The default for :data:`values` is to use the same list as
      specified for :data:`choices`.
    :param function OnChange: an optional routine that will be called
      when the  
      :class:`~wxpython.ComboBox` can be specified.
    :param (other): additional keyword arguments accepted by
      :class:`~wxpython.ComboBox` can be specified.
    '''
    def __init__(self,parent,dct,item,choices,values=None,OnChange=None,**kw):
        if values is None:
            values = choices
        if dct[item] in values:
            i = values.index(dct[item])
        else:
            i = 0
        startval = choices[i]
        wx.ComboBox.__init__(self,parent,wx.ID_ANY,startval,
                             choices = choices,
                             style=wx.CB_DROPDOWN|wx.CB_READONLY,
                             **kw)
        self.choices = choices
        self.values = values
        self.dct = dct
        self.item = item
        self.OnChange = OnChange
        self.Bind(wx.EVT_COMBOBOX, self.onSelection)
    def onSelection(self,event):
        # respond to a selection by setting the enum value in the CIF dictionary
        if self.GetValue() in self.choices: # should always be true!
            self.dct[self.item] = self.values[self.choices.index(self.GetValue())]
        else:
            self.dct[self.item] = self.values[0] # unknown
        if self.OnChange: self.OnChange(event)

################################################################################
class popupSelectorButton(wx.Button):
    '''Create a button that will invoke a menu with choices that can 
    be selected. Do special stuff if the first item is "all"

    TODO: It might be better to make this a wx.ComboCtrl if I can figure out
    how to make that work, or perhaps make that an option

    :param wx.Frame parent: a panel or frame that is the parent to this 
      button
    :param str lbl: a label for the button
    :param list choices: a list of str's with labels for the items in the
      menu
    :param list selected: a list of bool's that determine if the menu item
      is initial selected
    :param: dict choiceDict: a dict with both choices and their 
      values (selections). Use this or choices & selected, not both. 
      If this is used, the values are set as radiobutton choices, 
      only the most recent setting is selected. 
    :param function OnChange: an optional function that is called after the 
      menu is removed
    :param others: other keyword parameters are allowed. They will be 
      passed to the OnChange routine. 
    '''
    def __init__(self,parent,lbl,choices=None,selected=None,
                     choiceDict=None, OnChange=None,**kw):
        self.parent = parent
        self.choices = choices
        self.selected = selected
        self.choiceDict = None
        if choiceDict is not None:
            self.choices = list(choiceDict.keys())
            self.selected = list(choiceDict.values())
            self.choiceDict = choiceDict
        self.OnChange = OnChange
        self.kw = kw
        wx.Button.__init__(self,parent,label=lbl)
        self.Bind(wx.EVT_BUTTON, self.popupSelector)
    
    def popupSelector(self,event):
        '''Show the menu and then get current values. Optionally call the 
        OnChange routine. 
        '''
        menu = wx.Menu()
        menuList = []
        #breakpoint()
        for i in range(len(self.choices)):
            menuList.append(
                menu.Append(wx.ID_ANY, self.choices[i], '', wx.ITEM_CHECK)
                )
            if self.selected[i]: menuList[-1].Check()
        self.parent.PopupMenu(menu)
        new = -1
        for i,m in enumerate(menuList):
            if self.selected[i] != m.IsChecked(): new = i
            self.selected[i] = m.IsChecked()
        if self.choiceDict is not None:
            self.selected = len(self.selected) * [False]
            for i in self.choiceDict: self.choiceDict[i] = False
            self.choiceDict[self.choices[new]] = True
            self.selected[new] = True
        elif self.choices[0] == "all":
            # special handling when 1st item is all: turn on everything
            # when all is selected.
            if new == 0:
                self.selected[:] = [False] + (len(self.selected)-1)*[True]
                
        menu.Destroy()
        if self.OnChange: wx.CallAfter(self.OnChange,**self.kw)

################################################################################
class G2ChoiceButton(wx.Choice):
    '''A customized version of a wx.Choice that automatically initializes
    the control to match a supplied value and saves the choice directly
    into an array or list. Optionally a function can be called each time a
    choice is selected. The widget can be used with an array item that is set to 
    to the choice by number (``indLoc[indKey]``) or by string value
    (``strLoc[strKey]``) or both. The initial value is taken from ``indLoc[indKey]``
    if not None or ``strLoc[strKey]`` if not None. 

    :param wx.Panel parent: name of panel or frame that will be
      the parent to the widget. Can be None.
    :param list choiceList: a list or tuple of choices to offer the user.
    :param dict/list indLoc: a dict or list with the initial value to be
      placed in the Choice button. If this is None, this is ignored. 
    :param int/str indKey: the dict key or the list index for the value to be
      edited by the Choice button. If indLoc is not None then this
      must be specified and the ``indLoc[indKey]`` will be set. If the value
      for ``indLoc[indKey]`` is not None, it should be an integer in
      range(len(choiceList)). The Choice button will be initialized to the
      choice corresponding to the value in this element if not None.
    :param dict/list strLoc: a dict or list with the string value corresponding to
      indLoc/indKey. Default (None) means that this is not used. 
    :param int/str strKey: the dict key or the list index for the string value 
      The ``strLoc[strKey]`` element must exist or strLoc must be None (default).
    :param function onChoice: name of a function to call when the choice is made.
    '''
    def __init__(self,parent,choiceList,indLoc=None,indKey=None,strLoc=None,strKey=None,
                 onChoice=None,**kwargs):
        wx.Choice.__init__(self,parent,choices=choiceList,id=wx.ID_ANY,**kwargs)
        self.choiceList = choiceList
        self.indLoc = indLoc
        self.indKey = indKey
        self.strLoc = strLoc
        self.strKey = strKey
        self.onChoice = None
        self.SetSelection(wx.NOT_FOUND)
        if self.indLoc is not None and self.indKey is not None:
            try:
                self.SetSelection(self.indLoc[self.indKey])
                if self.strLoc is not None and self.strKey is not None:
                    self.strLoc[self.strKey] = self.GetStringSelection()
            except (KeyError,ValueError,TypeError):
                pass
        elif self.strLoc is not None and self.strKey is not None:
            try:
                self.SetSelection(choiceList.index(self.strLoc[self.strKey]))
                if self.indLoc is not None:
                    self.indLoc[self.indKey] = self.GetSelection()
            except (KeyError,ValueError,TypeError):
                pass
        self.Bind(wx.EVT_CHOICE, self._OnChoice)
        #if self.strLoc is not None: # make sure strLoc gets initialized
        #    self._OnChoice(None) # note that onChoice will not be called
        self.onChoice = onChoice
    def _OnChoice(self,event):
        if self.indLoc is not None:
            self.indLoc[self.indKey] = self.GetSelection()
        if self.strLoc is not None:
            self.strLoc[self.strKey] = self.GetStringSelection()
        if self.onChoice:
            self.onChoice()
    def setByString(self,string):
        'Find an entry matching string and select it'
        num = self.FindString(string)
        if num >= 0: self.SetSelection(num)

#### Custom checkbox that saves values into dict/list as used ##############################################################
class G2CheckBox(wx.CheckBox):
    '''A customized version of a CheckBox that automatically initializes
    the control to a supplied list or dict entry and updates that
    entry as the widget is used.

    :param wx.Panel parent: name of panel or frame that will be
      the parent to the widget. Can be None.
    :param str label: text to put on check button
    :param dict/list loc: the dict or list with the initial value to be
      placed in the CheckBox. 
    :param int/str key: the dict key or the list index for the value to be
      edited by the CheckBox. The ``loc[key]`` element must exist.
      The CheckBox will be initialized from this value.
      If the value is anything other that True (or 1), it will be taken as
      False.
    :param function OnChange: specifies a function or method that will be
      called when the CheckBox is changed (Default is None). 
      The called function is supplied with one argument, the calling event.
    '''
    def __init__(self,parent,label,loc,key,OnChange=None):
        wx.CheckBox.__init__(self,parent,id=wx.ID_ANY,label=label)
        self.loc = loc
        self.key = key
        self.OnChange = OnChange
        self.SetValue(self.loc[self.key]==True)
        self.Bind(wx.EVT_CHECKBOX, self._OnCheckBox)
    def _OnCheckBox(self,event):
        self.loc[self.key] = self.GetValue()
        if self.OnChange: self.OnChange(event)

def G2CheckBoxFrontLbl(parent,label,loc,key,OnChange=None):
    '''A customized version of a CheckBox that automatically initializes
    the control to a supplied list or dict entry and updates that
    entry as the widget is used. Same as :class:`G2CheckBox` except the 
    label is placed before the CheckBox and returns a sizer rather than the
    G2CheckBox. 

    If the CheckBox is needed, reference Sizer.myCheckBox.
    '''
    Sizer = wx.BoxSizer(wx.HORIZONTAL)
    Sizer.Add(wx.StaticText(parent,label=label),0,WACV)
    checkBox = G2CheckBox(parent,'',loc,key,OnChange)
    Sizer.Add(checkBox,0,WACV)
    Sizer.myCheckBox = checkBox
    return Sizer
    
def G2RadioButtons(parent,loc,key,choices,values=None,OnChange=None):
    '''A customized version of wx.RadioButton that returns a list 
    of coupled RadioButtons

    :param wx.Panel parent: name of panel or frame that will be
      the parent to the widgets. Can be None.
    :param dict/list loc: the dict or list with the initial value to be
      placed in the CheckBox. 
    :param int/str key: the dict key or the list index for the value to be
      edited by the CheckBox. The ``loc[key]`` element must exist.
      The CheckButton will be initialized from this value.
    :param list choices:
    :param list values:
    :param function OnChange: specifies a function or method that will be
      called when the CheckBox is changed (Default is None). 
      The called function is supplied with one argument, the calling event.
    '''
    def _OnEvent(event):
        if event.GetEventObject() not in buttons:
            if GSASIIpath.GetConfigValue('debug'): print('Strange: unknown button')
            return
        loc[key] = values[buttons.index(event.GetEventObject())]
        if OnChange: OnChange(event)
    if not values:
        values = list(range(len(choices)))
    buttons = []
    kw = {'style':wx.RB_GROUP}
    for i,c in enumerate(choices):
        if i == 1:
            kw = {}
        buttons.append(wx.RadioButton(parent,wx.ID_ANY,c,**kw))
        if loc[key] == values[i]: buttons[-1].SetValue(True)
        buttons[-1].Bind(wx.EVT_RADIOBUTTON, _OnEvent)
    return buttons

#### Commonly used dialogs ################################################################################
def CallScrolledMultiEditor(parent,dictlst,elemlst,prelbl=[],postlbl=[],
                 title='Edit items',header='',size=(300,250),
                             CopyButton=False, ASCIIonly=False, **kw):
    '''Shell routine to call a ScrolledMultiEditor dialog. See
    :class:`ScrolledMultiEditor` for parameter definitions.

    :returns: True if the OK button is pressed; False if the window is closed
      with the system menu or the Cancel button.

    '''
    dlg = ScrolledMultiEditor(parent,dictlst,elemlst,prelbl,postlbl,
                              title,header,size,
                              CopyButton, ASCIIonly, **kw)
    dlg.CenterOnParent()
    if dlg.ShowModal() == wx.ID_OK:
        dlg.Destroy()
        return True
    else:
        dlg.Destroy()
        return False

################################################################################
class ScrolledMultiEditor(wx.Dialog):
    '''Define a window for editing a potentially large number of dict- or
    list-contained values with validation for each item. Edited values are
    automatically placed in their source location. If invalid entries
    are provided, the TextCtrl is turned yellow and the OK button is disabled.

    The type for each TextCtrl validation is determined by the
    initial value of the entry (int, float or string). 
    Float values can be entered in the TextCtrl as numbers or also
    as algebraic expressions using operators + - / \\* () and \\*\\*,
    in addition pi, sind(), cosd(), tand(), and sqrt() can be used,
    as well as appreviations s(), sin(), c(), cos(), t(), tan() and sq(). 

    :param wx.Frame parent: name of parent window, or may be None

    :param tuple dictlst: a list of dicts or lists containing values to edit

    :param tuple elemlst: a list of keys/indices for items in dictlst. 
      Note that elemlst must have the same length as dictlst, where each 
      item in elemlst will will match an entry for an entry for successive
      dicts/lists in dictlst.

    :param tuple prelbl: a list of labels placed before the TextCtrl for each
      item (optional)
   
    :param tuple postlbl: a list of labels placed after the TextCtrl for each
      item (optional)

    :param str title: a title to place in the frame of the dialog

    :param str header: text to place at the top of the window. May contain
      new line characters. 

    :param wx.Size size: a size parameter that dictates the
      size for the scrolled region of the dialog. The default is
      (300,250). 

    :param bool CopyButton: if True adds a small button that copies the
      value for the current row to all fields below (default is False)
      
    :param bool ASCIIonly: if set as True will remove unicode characters from
      strings
      
    :param list minvals: optional list of minimum values for validation
      of float or int values. Ignored if value is None.
    :param list maxvals: optional list of maximum values for validation
      of float or int values. Ignored if value is None.
    :param list sizevals: optional list of wx.Size values for each input
      widget. Ignored if value is None.
      
    :param tuple checkdictlst: an optional list of dicts or lists containing bool
      values (similar to dictlst). 
    :param tuple checkelemlst: an optional list of dicts or lists containing bool
      key values (similar to elemlst). Must be used with checkdictlst.
    :param string checklabel: a string to use for each checkbutton
      
    :returns: the wx.Dialog created here. Use method .ShowModal() to display it.
    
    *Example for use of ScrolledMultiEditor:*

    ::

        dlg = <pkg>.ScrolledMultiEditor(frame,dictlst,elemlst,prelbl,postlbl,
                                        header=header)
        if dlg.ShowModal() == wx.ID_OK:
             for d,k in zip(dictlst,elemlst):
                 print d[k]

    *Example definitions for dictlst and elemlst:*

    ::
      
          dictlst = (dict1,list1,dict1,list1)
          elemlst = ('a', 1, 2, 3)

      This causes items dict1['a'], list1[1], dict1[2] and list1[3] to be edited.
    
    Note that these items must have int, float or str values assigned to
    them. The dialog will force these types to be retained. String values
    that are blank are marked as invalid. 
    '''
    
    def __init__(self,parent,dictlst,elemlst,prelbl=[],postlbl=[],
                 title='Edit items',header='',size=(300,250),
                 CopyButton=False, ASCIIonly=False,
                 minvals=[],maxvals=[],sizevals=[],
                 checkdictlst=[], checkelemlst=[], checklabel=""):
        if len(dictlst) != len(elemlst):
            raise Exception("ScrolledMultiEditor error: len(dictlst) != len(elemlst) "+str(len(dictlst))+" != "+str(len(elemlst)))
        if len(checkdictlst) != len(checkelemlst):
            raise Exception("ScrolledMultiEditor error: len(checkdictlst) != len(checkelemlst) "+str(len(checkdictlst))+" != "+str(len(checkelemlst)))
        wx.Dialog.__init__( # create dialog & sizer
            self,parent,wx.ID_ANY,title,
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.orig = []
        self.dictlst = dictlst
        self.elemlst = elemlst
        self.checkdictlst = checkdictlst
        self.checkelemlst = checkelemlst
        self.StartCheckValues = [checkdictlst[i][checkelemlst[i]] for i in range(len(checkdictlst))]
        self.ButtonIndex = {}
        for d,i in zip(dictlst,elemlst):
            self.orig.append(d[i])
        # add a header if supplied
        if header:
            subSizer = wx.BoxSizer(wx.HORIZONTAL)
            subSizer.Add((-1,-1),1,wx.EXPAND)
            subSizer.Add(wx.StaticText(self,wx.ID_ANY,header))
            subSizer.Add((-1,-1),1,wx.EXPAND)
            mainSizer.Add(subSizer,0,wx.EXPAND,0)
        # make OK button now, because we will need it for validation
        self.OKbtn = wx.Button(self, wx.ID_OK)
        self.OKbtn.SetDefault()
        # create scrolled panel and sizer
        panel = wxscroll.ScrolledPanel(self, wx.ID_ANY,size=size,
            style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        cols = 4
        if CopyButton: cols += 1
        subSizer = wx.FlexGridSizer(cols=cols,hgap=2,vgap=2)
        self.ValidatedControlsList = [] # make list of TextCtrls
        self.CheckControlsList = [] # make list of CheckBoxes
        for i,(d,k) in enumerate(zip(dictlst,elemlst)):
            if i >= len(prelbl): # label before TextCtrl, or put in a blank
                subSizer.Add((-1,-1)) 
            else:
                subSizer.Add(wx.StaticText(panel,wx.ID_ANY,str(prelbl[i])))
            kargs = {}
            if i < len(minvals):
                if minvals[i] is not None: kargs['xmin']=minvals[i]
            if i < len(maxvals):
                if maxvals[i] is not None: kargs['xmax']=maxvals[i]
            if i < len(sizevals):
                if sizevals[i]: kargs['size']=sizevals[i]
            if CopyButton:
                if i+1 == len(dictlst):
                    but = (-1,-1)
                else:
                    import wx.lib.colourselect as wscs  # is there a way to test? 
                    but = wscs.ColourSelect(label='v', # would like to use u'\u2193' or u'\u25BC' but not in WinXP
                                            parent=panel,colour=(255,255,200),size=wx.Size(30,23),
                                            style=wx.RAISED_BORDER)
                    but.Bind(wx.EVT_BUTTON, self._OnCopyButton)
                    if 'phoenix' in wx.version():
                        but.SetToolTip('Press to copy adjacent value to all rows below')
                    else:
                        but.SetToolTipString('Press to copy adjacent value to all rows below')
                    self.ButtonIndex[but] = i
                subSizer.Add(but)
            # create the validated TextCrtl, store it and add it to the sizer
            ctrl = ValidatedTxtCtrl(panel,d,k,OKcontrol=self.ControlOKButton,ASCIIonly=ASCIIonly,
                                    **kargs)
            self.ValidatedControlsList.append(ctrl)
            subSizer.Add(ctrl)
            if i < len(postlbl): # label after TextCtrl, or put in a blank
                subSizer.Add(wx.StaticText(panel,wx.ID_ANY,str(postlbl[i])))
            else:
                subSizer.Add((-1,-1))
            if i < len(checkdictlst):
                ch = G2CheckBox(panel,checklabel,checkdictlst[i],checkelemlst[i])
                self.CheckControlsList.append(ch)
                subSizer.Add(ch)                    
            else:
                subSizer.Add((-1,-1))
        # finish up ScrolledPanel
        panel.SetSizer(subSizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
        # patch for wx 2.9 on Mac
        i,j= wx.__version__.split('.')[0:2]
        if int(i)+int(j)/10. > 2.8 and 'wxOSX' in wx.PlatformInfo:
            panel.SetMinSize((subSizer.GetSize()[0]+30,panel.GetSize()[1]))        
        mainSizer.Add(panel,1, wx.ALL|wx.EXPAND,1)

        # Sizer for OK/Close buttons. N.B. on Close changes are discarded
        # by restoring the initial values
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btnsizer.Add(self.OKbtn)
        self.OKbtn.Bind(wx.EVT_BUTTON,lambda event: self.EndModal(wx.ID_OK))
        btn = wx.Button(self, wx.ID_CLOSE,"Cancel") 
        btn.Bind(wx.EVT_BUTTON,self._onClose)
        btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        # size out the window. Set it to be enlarged but not made smaller
        self.SetSizer(mainSizer)
        mainSizer.Fit(self)
        self.SetMinSize(self.GetSize())

    def _OnCopyButton(self,event):
        'Implements the copy down functionality'
        but = event.GetEventObject()
        n = self.ButtonIndex.get(but)
        if n is None: return
        for i,(d,k,ctrl) in enumerate(zip(self.dictlst,self.elemlst,self.ValidatedControlsList)):
            if i < n: continue
            if i == n:
                val = d[k]
                continue
            d[k] = val
            ctrl.SetValue(val)
        for i in range(len(self.checkdictlst)):
            if i < n: continue
            self.checkdictlst[i][self.checkelemlst[i]] = self.checkdictlst[n][self.checkelemlst[n]]
            self.CheckControlsList[i].SetValue(self.checkdictlst[i][self.checkelemlst[i]])
    def _onClose(self,event):
        'Used on Cancel: Restore original values & close the window'
        for d,i,v in zip(self.dictlst,self.elemlst,self.orig):
            d[i] = v
        for i in range(len(self.checkdictlst)):
            self.checkdictlst[i][self.checkelemlst[i]] = self.StartCheckValues[i]
        self.EndModal(wx.ID_CANCEL)
        
    def ControlOKButton(self,setvalue):
        '''Enable or Disable the OK button for the dialog. Note that this is
        passed into the ValidatedTxtCtrl for use by validators.

        :param bool setvalue: if True, all entries in the dialog are
          checked for validity. if False then the OK button is disabled.

        '''
        if setvalue: # turn button on, do only if all controls show as valid
            for ctrl in self.ValidatedControlsList:
                if ctrl.invalid:
                    self.OKbtn.Disable()
                    return
            else:
                self.OKbtn.Enable()
        else:
            self.OKbtn.Disable()

###############################################  Multichoice Dialog with set all, toggle & filter options
class G2MultiChoiceDialog(wx.Dialog):
    '''A dialog similar to wx.MultiChoiceDialog except that buttons are
    added to set all choices and to toggle all choices and a filter is 
    available to select from available entries. Note that if multiple
    entries are placed in the filter box separated by spaces, all 
    of the strings must be present for an item to be shown. 

    :param wx.Frame ParentFrame: reference to parent frame
    :param str title: heading above list of choices
    :param str header: Title to place on window frame 
    :param list ChoiceList: a list of choices where one more will be selected
    :param bool toggle: If True (default) the toggle and select all buttons
      are displayed
    :param bool monoFont: If False (default), use a variable-spaced font;
      if True use a equally-spaced font.
    :param bool filterBox: If True (default) an input widget is placed on
      the window and only entries matching the entered text are shown.
    :param dict extraOpts: a dict containing a entries of form label_i and value_i with extra
      options to present to the user, where value_i is the default value. 
      Options are listed ordered by the value_i values.
    :param list selected: list of indicies for items that should be 
    :param kw: optional keyword parameters for the wx.Dialog may
      be included such as size [which defaults to `(320,310)`] and
      style (which defaults to `wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.CENTRE| wx.OK | wx.CANCEL`);
      note that `wx.OK` and `wx.CANCEL` style items control 
      the presence of the eponymous buttons in the dialog.
    :returns: the name of the created dialog  
    '''
    def __init__(self,parent, title, header, ChoiceList, toggle=True,
                 monoFont=False, filterBox=True, extraOpts={}, selected=[],
                 **kw):
        # process keyword parameters, notably style
        options = {'size':(320,310), # default Frame keywords
                   'style':wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.CENTRE| wx.OK | wx.CANCEL,
                   }
        options.update(kw)
        self.ChoiceList = ['%4d) %s'%(i,item) for i,item in enumerate(ChoiceList)] # numbered list of choices (list of str values)
        self.Selections = len(self.ChoiceList) * [False,] # selection status for each choice (list of bools)
        for i in selected:
            self.Selections[i] = True
        self.filterlist = range(len(self.ChoiceList)) # list of the choice numbers that have been filtered (list of int indices)
        self.Stride = 1
        if options['style'] & wx.OK:
            useOK = True
            options['style'] ^= wx.OK
        else:
            useOK = False
        if options['style'] & wx.CANCEL:
            useCANCEL = True
            options['style'] ^= wx.CANCEL
        else:
            useCANCEL = False        
        # create the dialog frame
        wx.Dialog.__init__(self,parent,wx.ID_ANY,header,**options)
        # fill the dialog
        Sizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(self,wx.ID_ANY,title,size=(-1,35)),
            1,wx.ALL|wx.EXPAND,1)
        if filterBox:
            self.timer = wx.Timer()
            self.timer.Bind(wx.EVT_TIMER,self.Filter)
            topSizer.Add(wx.StaticText(self,wx.ID_ANY,'Name \nFilter: '),0,wx.ALL|WACV,1)
            self.filterBox = wx.TextCtrl(self, wx.ID_ANY, size=(80,-1),style=wx.TE_PROCESS_ENTER)
            self.filterBox.Bind(wx.EVT_TEXT,self.onChar)
            self.filterBox.Bind(wx.EVT_TEXT_ENTER,self.Filter)
            topSizer.Add(self.filterBox,0,wx.ALL|WACV,0)
        Sizer.Add(topSizer,0,wx.ALL|wx.EXPAND,8)
        self.settingRange = False
        self.rangeFirst = None
        self.clb = wx.CheckListBox(self, wx.ID_ANY, (30,30), wx.DefaultSize, self.ChoiceList)
        self._ShowSelections()
        self.clb.Bind(wx.EVT_CHECKLISTBOX,self.OnCheck)
        if monoFont:
            font1 = wx.Font(self.clb.GetFont().GetPointSize(),
                            wx.MODERN, wx.NORMAL, wx.NORMAL, False)
            self.clb.SetFont(font1)
        Sizer.Add(self.clb,1,wx.LEFT|wx.RIGHT|wx.EXPAND,10)
        Sizer.Add((-1,10))
        # set/toggle buttons
        if toggle:
            tSizer = wx.FlexGridSizer(cols=2,hgap=5,vgap=5)
            tSizer.Add(wx.StaticText(self,label=' Apply stride:'),0,WACV)
            numbs = [str(i+1) for i in range(9)]+[str(2*i+10) for i in range(6)]
            self.stride = wx.ComboBox(self,value='1',choices=numbs,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            self.stride.Bind(wx.EVT_COMBOBOX,self.OnStride)
            tSizer.Add(self.stride,0,WACV)
            setBut = wx.Button(self,wx.ID_ANY,'Set All')
            setBut.Bind(wx.EVT_BUTTON,self._SetAll)
            tSizer.Add(setBut)
            togBut = wx.Button(self,wx.ID_ANY,'Toggle All')
            togBut.Bind(wx.EVT_BUTTON,self._ToggleAll)
            tSizer.Add(togBut)
            self.rangeBut = wx.ToggleButton(self,wx.ID_ANY,'Set Range')
            self.rangeBut.Bind(wx.EVT_TOGGLEBUTTON,self.SetRange)
            tSizer.Add(self.rangeBut)
            self.rangeCapt = wx.StaticText(self,wx.ID_ANY,'')
            tSizer.Add(self.rangeCapt)
            Sizer.Add(tSizer,0,wx.LEFT,12)
        # Extra widgets
        Sizer.Add((-1,5),0,wx.LEFT,0)
        bSizer = wx.BoxSizer(wx.VERTICAL)
        for lbl in sorted(extraOpts.keys()):
            if not lbl.startswith('label'): continue
            key = lbl.replace('label','value')
            if key not in extraOpts: continue
            eSizer = wx.BoxSizer(wx.HORIZONTAL)
            if type(extraOpts[key]) is bool:
                eSizer.Add(G2CheckBox(self,extraOpts[lbl],extraOpts,key))
            else:
                eSizer.Add(wx.StaticText(self,wx.ID_ANY,extraOpts[lbl]))
                eSizer.Add(ValidatedTxtCtrl(self,extraOpts,key))
            bSizer.Add(eSizer,0,wx.LEFT,0)
        Sizer.Add(bSizer,0,wx.CENTER,0)
        Sizer.Add((-1,5),0,wx.LEFT,0)
        # OK/Cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        if useOK:
            self.OKbtn = wx.Button(self, wx.ID_OK)
            self.OKbtn.SetDefault()
            btnsizer.AddButton(self.OKbtn)
            self.OKbtn.Bind(wx.EVT_BUTTON,self.onOk)
        if useCANCEL:
            btn = wx.Button(self, wx.ID_CANCEL)
            btn.Bind(wx.EVT_BUTTON,self.onCancel)
            btnsizer.AddButton(btn)
        btnsizer.Realize()
        Sizer.Add((-1,5))
        Sizer.Add(btnsizer,0,wx.ALIGN_RIGHT,50)
        Sizer.Add((-1,20))
        # OK done, let's get outa here
        self.SetSizer(Sizer)
        Sizer.Fit(self)
        self.CenterOnParent()
        
    def onOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def onCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_CANCEL)
        
    def OnStride(self,event):
        self.Stride = int(self.stride.GetValue())

    def SetRange(self,event):
        '''Respond to a press of the Set Range button. Set the range flag and
        the caption next to the button
        '''
        self.settingRange = self.rangeBut.GetValue()
        if self.settingRange:
            self.rangeCapt.SetLabel('Select range start')
        else:
            self.rangeCapt.SetLabel('')            
        self.rangeFirst = None
        
    def GetSelections(self):
        'Returns a list of the indices for the selected choices'
        # update self.Selections with settings for displayed items
        for i in range(len(self.filterlist)):
            self.Selections[self.filterlist[i]] = self.clb.IsChecked(i)
        # return all selections, shown or hidden
        return [i for i in range(len(self.Selections)) if self.Selections[i]]
        
    def SetSelections(self,selList):
        '''Sets the selection indices in selList as selected. Resets any previous
        selections for compatibility with wx.MultiChoiceDialog. Note that
        the state for only the filtered items is shown.

        :param list selList: indices of items to be selected. These indices
          are referenced to the order in self.ChoiceList
        '''
        self.Selections = len(self.ChoiceList) * [False,] # reset selections
        for sel in selList:
            self.Selections[sel] = True
        self._ShowSelections()

    def _ShowSelections(self):
        'Show the selection state for displayed items'
        if 'phoenix' in wx.version():
            self.clb.SetCheckedItems(
                [i for i in range(len(self.filterlist)) if self.Selections[self.filterlist[i]]]
            ) # Note anything previously checked will be cleared.
        else:
            self.clb.SetChecked(
                [i for i in range(len(self.filterlist)) if self.Selections[self.filterlist[i]]]
            ) # Note anything previously checked will be cleared.
            
    def _SetAll(self,event):
        'Set all viewed choices on'
        if 'phoenix' in wx.version():
            self.clb.SetCheckedItems(range(0,len(self.filterlist),self.Stride))
        else:
            self.clb.SetChecked(range(0,len(self.filterlist),self.Stride))
        self.stride.SetValue('1')
        self.Stride = 1
        
    def _ToggleAll(self,event):
        'flip the state of all viewed choices'
        for i in range(len(self.filterlist)):
            self.clb.Check(i,not self.clb.IsChecked(i))
            
    def onChar(self,event):
        'Respond to keyboard events in the Filter box'
        self.OKbtn.Enable(False)
        if self.timer.IsRunning():
            self.timer.Stop()
        self.timer.Start(1000,oneShot=True)
        if event: event.Skip()
        
    def OnCheck(self,event):
        '''for CheckListBox events; if Set Range is in use, this sets/clears all
        entries in range between start and end according to the value in start.
        Repeated clicks on the start change the checkbox state, but do not trigger
        the range copy. 
        The caption next to the button is updated on the first button press.
        '''
        if self.settingRange:
            id = event.GetInt()
            if self.rangeFirst is None:
                name = self.clb.GetString(id)
                self.rangeCapt.SetLabel(name+' to...')
                self.rangeFirst = id
            elif self.rangeFirst == id:
                pass
            else:
                for i in range(min(self.rangeFirst,id), max(self.rangeFirst,id)+1,self.Stride):
                    self.clb.Check(i,self.clb.IsChecked(self.rangeFirst))
                self.rangeBut.SetValue(False)
                self.rangeCapt.SetLabel('')
            return
        
    def Filter(self,event):
        '''Read text from filter control and select entries that match. Called by
        Timer after a delay with no input or if Enter is pressed.
        '''
        if self.timer.IsRunning():
            self.timer.Stop()
        self.GetSelections() # record current selections
        txt = self.filterBox.GetValue()
        txt = txt.lower()
        self.clb.Clear()
        
        self.Update()
        self.filterlist = []
        if txt:
            ChoiceList = []
            for i,item in enumerate(self.ChoiceList):
                for t in txt.split():
                    if item.lower().find(t) == -1:
                        break
                else:
                    ChoiceList.append(item)
                    self.filterlist.append(i)
        else:
            self.filterlist = range(len(self.ChoiceList))
            ChoiceList = self.ChoiceList
        self.clb.AppendItems(ChoiceList)
        self._ShowSelections()
        self.OKbtn.Enable(True)

###############################################  Multichoice in a sizer with set all, toggle & filter options
class G2MultiChoiceWindow(wx.BoxSizer):
    '''Creates a sizer similar to G2MultiChoiceDialog except that 
    buttons are added to set all choices and to toggle all choices. This
    is placed in a sizer, so that it can be used in a frame or panel. 

    :param parent: reference to parent frame/panel
    :param str title: heading above list of choices
    :param list ChoiceList: a list of choices where one more will be selected
    :param list SelectList: a list of selected choices
    :param bool toggle: If True (default) the toggle and select all buttons
      are displayed
    :param bool monoFont: If False (default), use a variable-spaced font;
      if True use a equally-spaced font.
    :param bool filterBox: If True (default) an input widget is placed on
      the window and only entries matching the entered text are shown.
    :param function OnChange: a reference to a callable object, that 
      is called each time any a choice is changed. Default is None which
      will not be called. 
    :param list OnChangeArgs: a list of arguments to be supplied to function
      OnChange. The default is a null list.
    :returns: the name of the created sizer
    '''
    def __init__(self, parent, title, ChoiceList, SelectList, toggle=True,
                 monoFont=False, filterBox=True,
                     OnChange=None, OnChangeArgs=[], helpText=None):
        self.SelectList = SelectList
        self.ChoiceList = ['%4d) %s'%(i,item) for i,item in enumerate(ChoiceList)] # numbered list of choices (list of str values)
        self.frm = parent
        self.Selections = len(self.ChoiceList) * [False,] # selection status for each choice (list of bools)
        self.filterlist = range(len(self.ChoiceList)) # list of the choice numbers that have been filtered (list of int indices)
        self.Stride = 1
        self.OnChange = OnChange
        self.OnChangeArgs = OnChangeArgs
        # fill frame
        wx.BoxSizer.__init__(self,wx.VERTICAL)
        # fill the sizer
        Sizer = self
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(wx.StaticText(self.frm,wx.ID_ANY,title,size=(-1,35)),0,WACV)
        if helpText:
            topSizer.Add(HelpButton(self.frm,helpText,wrap=400),0,wx.ALL,5)
        topSizer.Add((1,-1),1,wx.ALL|wx.EXPAND,1)
        if filterBox:
            self.timer = wx.Timer()
            self.timer.Bind(wx.EVT_TIMER,self.Filter)
            topSizer.Add(wx.StaticText(self.frm,wx.ID_ANY,'Name \nFilter: '),0,wx.ALL|WACV,1)
            self.filterBox = wx.TextCtrl(self.frm, wx.ID_ANY, size=(80,-1),style=wx.TE_PROCESS_ENTER)
            self.filterBox.Bind(wx.EVT_TEXT,self.onChar)
            self.filterBox.Bind(wx.EVT_TEXT_ENTER,self.Filter)
            topSizer.Add(self.filterBox,0,wx.ALL|WACV,0)
        Sizer.Add(topSizer,0,wx.ALL|wx.EXPAND,8)
        self.settingRange = False
        self.rangeFirst = None
        self.clb = wx.CheckListBox(self.frm, wx.ID_ANY, (30,30), wx.DefaultSize, self.ChoiceList)
        self.clb.Bind(wx.EVT_CHECKLISTBOX,self.OnCheck)
        if monoFont:
            font1 = wx.Font(self.clb.GetFont().GetPointSize(),
                            wx.MODERN, wx.NORMAL, wx.NORMAL, False)
            self.clb.SetFont(font1)
        Sizer.Add(self.clb,1,wx.LEFT|wx.RIGHT|wx.EXPAND,10)
        Sizer.Add((-1,10))
        # set/toggle buttons
        if toggle:
            tSizer = wx.BoxSizer(wx.HORIZONTAL)
            tSizer.Add(wx.StaticText(self.frm,label=' Apply stride:'),0,WACV)
            numbs = [str(i+1) for i in range(9)]+[str(2*i+10) for i in range(6)]
            self.stride = wx.ComboBox(self.frm,value='1',choices=numbs,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            self.stride.Bind(wx.EVT_COMBOBOX,self.OnStride)
            tSizer.Add(self.stride,0,WACV)
            Sizer.Add(tSizer,0,wx.LEFT,12)
            tSizer = wx.BoxSizer(wx.HORIZONTAL)
            setBut = wx.Button(self.frm,wx.ID_ANY,'Set All')
            setBut.Bind(wx.EVT_BUTTON,self._SetAll)
            tSizer.Add(setBut)
            togBut = wx.Button(self.frm,wx.ID_ANY,'Toggle All')
            togBut.Bind(wx.EVT_BUTTON,self._ToggleAll)
            tSizer.Add(togBut)
            self.rangeBut = wx.ToggleButton(self.frm,wx.ID_ANY,'Set Range')
            self.rangeBut.Bind(wx.EVT_TOGGLEBUTTON,self.SetRange)
            tSizer.Add(self.rangeBut)
            Sizer.Add(tSizer,0,wx.LEFT,12)
            tSizer = wx.BoxSizer(wx.HORIZONTAL)
            self.rangeCapt = wx.StaticText(self.frm,wx.ID_ANY,'')
            tSizer.Add(self.rangeCapt,1,wx.EXPAND,1)
            Sizer.Add(tSizer,0,wx.LEFT,12)
        self.SetSelections(self.SelectList)
                
    def OnStride(self,event):
        self.Stride = int(self.stride.GetValue())

    def SetRange(self,event):
        '''Respond to a press of the Set Range button. Set the range flag and
        the caption next to the button
        '''
        self.settingRange = self.rangeBut.GetValue()
        if self.settingRange:
            self.rangeCapt.SetLabel('Select range start')
        else:
            self.rangeCapt.SetLabel('')            
        self.rangeFirst = None
        
    def GetSelections(self):
        'Returns a list of the indices for the selected choices'
        # update self.Selections with settings for displayed items
        for i in range(len(self.filterlist)):
            self.Selections[self.filterlist[i]] = self.clb.IsChecked(i)
        # return all selections, shown or hidden
        return [i for i in range(len(self.Selections)) if self.Selections[i]]
        
    def SetSelections(self,selList):
        '''Sets the selection indices in selList as selected. Resets any previous
        selections for compatibility with wx.MultiChoiceDialog. Note that
        the state for only the filtered items is shown.

        :param list selList: indices of items to be selected. These indices
          are referenced to the order in self.ChoiceList
        '''
        self.Selections = len(self.ChoiceList) * [False,] # reset selections
        for sel in selList:
            self.Selections[sel] = True
        self._ShowSelections()

    def _ShowSelections(self):
        'Show the selection state for displayed items'
        if 'phoenix' in wx.version():
            self.clb.SetCheckedItems(
                [i for i in range(len(self.filterlist)) if self.Selections[self.filterlist[i]]]
            ) # Note anything previously checked will be cleared.
        else:
            self.clb.SetChecked(
                [i for i in range(len(self.filterlist)) if self.Selections[self.filterlist[i]]]
            ) # Note anything previously checked will be cleared.
        if self.OnChange:
            self.OnChange(self.GetSelections(),*self.OnChangeArgs)
        try:
            self.SelectList.clear()
        except:  # patch: clear not in EPD
            for i in reversed((range(len(self.SelectList)))):
                del self.SelectList[i]
        for i,val in enumerate(self.Selections):
            if val: self.SelectList.append(i)

    def _SetAll(self,event):
        'Set all viewed choices on'
        if 'phoenix' in wx.version():
            self.clb.SetCheckedItems(range(0,len(self.filterlist),self.Stride))
        else:
            self.clb.SetChecked(range(0,len(self.filterlist),self.Stride))
        self.stride.SetValue('1')
        self.Stride = 1
        self.GetSelections() # record current selections
        self._ShowSelections()
        
    def _ToggleAll(self,event):
        'flip the state of all viewed choices'
        for i in range(len(self.filterlist)):
            self.clb.Check(i,not self.clb.IsChecked(i))
        self.GetSelections() # record current selections
        self._ShowSelections()
            
    def onChar(self,event):
        'Respond to keyboard events in the Filter box'
        if self.timer.IsRunning():
            self.timer.Stop()
        self.timer.Start(1000,oneShot=True)
        if event: event.Skip()
        
    def OnCheck(self,event):
        '''for CheckListBox events; if Set Range is in use, this sets/clears all
        entries in range between start and end according to the value in start.
        Repeated clicks on the start change the checkbox state, but do not trigger
        the range copy. 
        The caption next to the button is updated on the first button press.
        '''
        if self.settingRange:
            id = event.GetInt()
            if self.rangeFirst is None:
                name = self.clb.GetString(id)
                self.rangeCapt.SetLabel(name+' to...')
                self.rangeFirst = id
            elif self.rangeFirst == id:
                pass
            else:
                for i in range(min(self.rangeFirst,id), max(self.rangeFirst,id)+1,self.Stride):
                    self.clb.Check(i,self.clb.IsChecked(self.rangeFirst))
                self.rangeBut.SetValue(False)
                self.rangeCapt.SetLabel('')
                self.settingRange = False
                self.rangeFirst = None
        self.GetSelections() # record current selections
        self._ShowSelections()

    def Filter(self,event):
        '''Read text from filter control and select entries that match. Called by
        Timer after a delay with no input or if Enter is pressed.
        '''
        if self.timer.IsRunning():
            self.timer.Stop()
        self.GetSelections() # record current selections
        txt = self.filterBox.GetValue()
        self.clb.Clear()
        
        self.filterlist = []
        if txt:
            txt = txt.lower()
            ChoiceList = []
            for i,item in enumerate(self.ChoiceList):
                if item.lower().find(txt) != -1:
                    ChoiceList.append(item)
                    self.filterlist.append(i)
        else:
            self.filterlist = range(len(self.ChoiceList))
            ChoiceList = self.ChoiceList
        self.clb.AppendItems(ChoiceList)
        self._ShowSelections()
        
def SelectEdit1Var(G2frame,array,labelLst,elemKeysLst,dspLst,refFlgElem):
    '''Select a variable from a list, then edit it and select histograms
    to copy it to.

    :param wx.Frame G2frame: main GSAS-II frame
    :param dict array: the array (dict or list) where values to be edited are kept
    :param list labelLst: labels for each data item
    :param list elemKeysLst: a list of lists of keys needed to be applied (see below)
      to obtain the value of each parameter
    :param list dspLst: list list of digits to be displayed (10,4) is 10 digits
      with 4 decimal places. Can be None.
    :param list refFlgElem: a list of lists of keys needed to be applied (see below)
      to obtain the refine flag for each parameter or None if the parameter
      does not have refine flag.

    Example::
      array = data 
      labelLst = ['v1','v2']
      elemKeysLst = [['v1'], ['v2',0]]
      refFlgElem = [None, ['v2',1]]

     * The value for v1 will be in data['v1'] and this cannot be refined while,
     * The value for v2 will be in data['v2'][0] and its refinement flag is data['v2'][1]
    '''
    def unkey(dct,keylist):
        '''dive into a nested set of dicts/lists applying keys in keylist
        consecutively
        '''
        d = dct
        for k in keylist:
            d = d[k]
        return d

    def OnChoice(event):
        'Respond when a parameter is selected in the Choice box'
        if 'phoenix' in wx.version():
            valSizer.Clear(True)
        else:
            valSizer.DeleteWindows()
        lbl = event.GetString()
        copyopts['currentsel'] = lbl
        i = labelLst.index(lbl)
        OKbtn.Enable(True)
        ch.SetLabel(lbl)
        args = {}
        if dspLst[i]:
            args = {'nDig':dspLst[i]}
        Val = ValidatedTxtCtrl(
            dlg,
            unkey(array,elemKeysLst[i][:-1]),
            elemKeysLst[i][-1],
            **args)
        copyopts['startvalue'] = unkey(array,elemKeysLst[i])
        #unkey(array,elemKeysLst[i][:-1])[elemKeysLst[i][-1]] = 
        valSizer.Add(Val,0,wx.LEFT,5)
        dlg.SendSizeEvent()
        
    # SelectEdit1Var execution begins here
    saveArray = copy.deepcopy(array) # keep original values
    TreeItemType = G2frame.GPXtree.GetItemText(G2frame.PickId)
    copyopts = {'InTable':False,"startvalue":None,'currentsel':None}        
    hst = G2frame.GPXtree.GetItemText(G2frame.PatternId)
    histList = G2pdG.GetHistsLikeSelected(G2frame)
    if not histList:
        G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame)
        return
    dlg = wx.Dialog(G2frame,wx.ID_ANY,'Set a parameter value',
        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5))
    subSizer = wx.BoxSizer(wx.HORIZONTAL)
    subSizer.Add((-1,-1),1,wx.EXPAND)
    subSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Select a parameter and set a new value'))
    subSizer.Add((-1,-1),1,wx.EXPAND)
    mainSizer.Add(subSizer,0,wx.EXPAND,0)
    mainSizer.Add((0,10))

    subSizer = wx.FlexGridSizer(0,2,5,0)
    subSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Parameter: '))
    ch = wx.Choice(dlg, wx.ID_ANY, choices = sorted(labelLst))
    ch.SetSelection(-1)
    ch.Bind(wx.EVT_CHOICE, OnChoice)
    subSizer.Add(ch)
    subSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Value: '))
    valSizer = wx.BoxSizer(wx.HORIZONTAL)
    subSizer.Add(valSizer)
    mainSizer.Add(subSizer)

    mainSizer.Add((-1,20))
    subSizer = wx.BoxSizer(wx.HORIZONTAL)
    subSizer.Add(G2CheckBox(dlg, 'Edit in table ', copyopts, 'InTable'))
    mainSizer.Add(subSizer)

    btnsizer = wx.StdDialogButtonSizer()
    OKbtn = wx.Button(dlg, wx.ID_OK,'Continue')
    OKbtn.Enable(False)
    OKbtn.SetDefault()
    OKbtn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_OK))
    btnsizer.AddButton(OKbtn)
    btn = wx.Button(dlg, wx.ID_CANCEL)
    btnsizer.AddButton(btn)
    btnsizer.Realize()
    mainSizer.Add((-1,5),1,wx.EXPAND,1)
    mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER,0)
    mainSizer.Add((-1,10))

    dlg.SetSizer(mainSizer)
    dlg.CenterOnParent()
    if dlg.ShowModal() != wx.ID_OK:
        array.update(saveArray)
        dlg.Destroy()
        return
    dlg.Destroy()

    copyList = []
    lbl = copyopts['currentsel']
    dlg = G2MultiChoiceDialog(G2frame,'Copy parameter '+lbl+' from\n'+hst,
        'Copy parameters', histList)
    dlg.CenterOnParent()
    try:
        if dlg.ShowModal() == wx.ID_OK:
            for i in dlg.GetSelections(): 
                copyList.append(histList[i])
        else:
            # reset the parameter since cancel was pressed
            array.update(saveArray)
            return
    finally:
        dlg.Destroy()

    prelbl = [hst]
    i = labelLst.index(lbl)
    keyLst = elemKeysLst[i]
    refkeys = refFlgElem[i]
    dictlst = [unkey(array,keyLst[:-1])]
    if refkeys is not None:
        refdictlst = [unkey(array,refkeys[:-1])]
    else:
        refdictlst = None
    Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,hst)
    hstData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters'))[0]
    for h in copyList:
        Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,h)
        instData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Instrument Parameters'))[0]
        if len(hstData) != len(instData) or hstData['Type'][0] != instData['Type'][0]:  #don't mix data types or lam & lam1/lam2 parms!
            print (h+' not copied - instrument parameters not commensurate')
            continue
        hData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,TreeItemType))
        if TreeItemType == 'Instrument Parameters':
            hData = hData[0]
        #copy the value if it is changed or we will not edit in a table
        valNow = unkey(array,keyLst)
        if copyopts['startvalue'] != valNow or not copyopts['InTable']:
            unkey(hData,keyLst[:-1])[keyLst[-1]] = valNow
        prelbl += [h]
        dictlst += [unkey(hData,keyLst[:-1])]
        if refdictlst is not None:
            refdictlst += [unkey(hData,refkeys[:-1])]
    if refdictlst is None:
        args = {}
    else:
        args = {'checkdictlst':refdictlst,
                'checkelemlst':len(dictlst)*[refkeys[-1]],
                'checklabel':'Refine?'}
    if copyopts['InTable']:
        dlg = ScrolledMultiEditor(
            G2frame,dictlst,
            len(dictlst)*[keyLst[-1]],prelbl,
            header='Editing parameter '+lbl,
            CopyButton=True,**args)
        dlg.CenterOnParent()
        if dlg.ShowModal() != wx.ID_OK:
            array.update(saveArray)
        dlg.Destroy()

#####  Single choice Dialog with filter options ###############################################################
class G2SingleChoiceDialog(wx.Dialog):
    '''A dialog similar to wx.SingleChoiceDialog except that a filter can be
    added.

    :param wx.Frame ParentFrame: reference to parent frame
    :param str title: heading above list of choices
    :param str header: Title to place on window frame 
    :param list ChoiceList: a list of choices where one will be selected
    :param bool monoFont: If False (default), use a variable-spaced font;
      if True use a equally-spaced font.
    :param bool filterBox: If True (default) an input widget is placed on
      the window and only entries matching the entered text are shown.
    :param kw: optional keyword parameters for the wx.Dialog may
      be included such as size [which defaults to `(320,310)`] and
      style (which defaults to ``wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.CENTRE | wx.OK | wx.CANCEL``);
      note that ``wx.OK`` and ``wx.CANCEL`` controls
      the presence of the eponymous buttons in the dialog.
    :returns: the name of the created dialog

    Example::

            dlg = G2SingleChoiceDialog(G2frame,'Select option from list',
                                           'Select option',optList)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    sel = optList[dlg.GetSelection()]
                else:
                    return
            finally:
                dlg.Destroy()

    '''
    def __init__(self,parent, title, header, ChoiceList, 
                 monoFont=False, filterBox=True, **kw):
        # process keyword parameters, notably style
        options = {'size':(320,310), # default Frame keywords
                   'style':wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.CENTRE| wx.OK | wx.CANCEL,
                   }
        options.update(kw)
        self.ChoiceList = ChoiceList
        self.filterlist = range(len(self.ChoiceList))
        if options['style'] & wx.OK:
            useOK = True
            options['style'] ^= wx.OK
        else:
            useOK = False
        if options['style'] & wx.CANCEL:
            useCANCEL = True
            options['style'] ^= wx.CANCEL
        else:
            useCANCEL = False        
        # create the dialog frame
        wx.Dialog.__init__(self,parent,wx.ID_ANY,header,**options)
        # fill the dialog
        Sizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        h = max(35,17*int(len(title)/26.+1)) # adjust height of title box with guessed # of lines
        topSizer.Add(wx.StaticText(self,wx.ID_ANY,title,size=(-1,h)),
            1,wx.ALL|wx.EXPAND,1)
        if filterBox:
            self.timer = wx.Timer()
            self.timer.Bind(wx.EVT_TIMER,self.Filter)
            topSizer.Add(wx.StaticText(self,wx.ID_ANY,'Filter: '),0,wx.ALL,1)
            self.filterBox = wx.TextCtrl(self, wx.ID_ANY, size=(80,-1),style=wx.TE_PROCESS_ENTER)
            self.filterBox.Bind(wx.EVT_CHAR,self.onChar)
            self.filterBox.Bind(wx.EVT_TEXT_ENTER,self.Filter)
            topSizer.Add(self.filterBox,0,wx.ALL,0)
        Sizer.Add(topSizer,0,wx.ALL|wx.EXPAND,8)
        self.clb = wx.ListBox(self, wx.ID_ANY, (30,30), wx.DefaultSize, ChoiceList)
        self.clb.Bind(wx.EVT_LEFT_DCLICK,self.onDoubleClick)
        if monoFont:
            font1 = wx.Font(self.clb.GetFont().GetPointSize(),wx.MODERN, wx.NORMAL, wx.NORMAL, False)
            self.clb.SetFont(font1)
        Sizer.Add(self.clb,1,wx.LEFT|wx.RIGHT|wx.EXPAND,10)
        Sizer.Add((-1,10))
        # OK/Cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        if useOK:
            self.OKbtn = wx.Button(self, wx.ID_OK)
            self.OKbtn.SetDefault()
            btnsizer.AddButton(self.OKbtn)
        if useCANCEL:
            btn = wx.Button(self, wx.ID_CANCEL)
            btnsizer.AddButton(btn)
        btnsizer.Realize()
        Sizer.Add((-1,5))
        Sizer.Add(btnsizer,0,wx.ALIGN_RIGHT,50)
        Sizer.Add((-1,20))
        # OK done, let's get outa here
        self.SetSizer(Sizer)
    def GetSelection(self):
        'Returns the index of the selected choice'
        i = self.clb.GetSelection()
        if i < 0 or i >= len(self.filterlist):
            return wx.NOT_FOUND
        return self.filterlist[i]
    def onChar(self,event):
        self.OKbtn.Enable(False)
        if self.timer.IsRunning():
            self.timer.Stop()
        self.timer.Start(1000,oneShot=True)
        if event: event.Skip()
    def Filter(self,event):
        if self.timer.IsRunning():
            self.timer.Stop()
        txt = self.filterBox.GetValue()
        self.clb.Clear()
        self.Update()
        self.filterlist = []
        if txt:
            txt = txt.lower()
            ChoiceList = []
            for i,item in enumerate(self.ChoiceList):
                if item.lower().find(txt) != -1:
                    ChoiceList.append(item)
                    self.filterlist.append(i)
        else:
            self.filterlist = range(len(self.ChoiceList))
            ChoiceList = self.ChoiceList
        self.clb.AppendItems(ChoiceList)
        self.OKbtn.Enable(True)
    def onDoubleClick(self,event):
        self.EndModal(wx.ID_OK)
        
################################################################################
class FlagSetDialog(wx.Dialog):
    ''' Creates popup with table of variables to be checked for e.g. refinement flags
    '''
    def __init__(self,parent,title,colnames,rownames,flags):
        wx.Dialog.__init__(self,parent,-1,title,
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = None
        self.colnames = colnames
        self.rownames = rownames
        self.flags = flags
        self.newflags = copy.copy(flags)
        self.Draw()
        
    def Draw(self):
        Indx = {}
        
        def OnSelection(event):
            Obj = event.GetEventObject()
            [name,ia] = Indx[Obj.GetId()]
            self.newflags[name][ia] = Obj.GetValue()
            
        if self.panel:
            self.panel.DestroyChildren()  #safe: wx.Panel
            self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        flagSizer = wx.FlexGridSizer(0,len(self.colnames),5,5)
        for item in self.colnames:
            flagSizer.Add(wx.StaticText(self.panel,label=item),0,WACV)
        for ia,atm in enumerate(self.rownames):
            flagSizer.Add(wx.StaticText(self.panel,label=atm),0,WACV)
            for name in self.colnames[1:]:
                if self.flags[name][ia]:
                    self.newflags[name][ia] = False     #default is off
                    flg = wx.CheckBox(self.panel,-1,label='')
                    flg.Bind(wx.EVT_CHECKBOX,OnSelection)
                    Indx[flg.GetId()] = [name,ia]
                    flagSizer.Add(flg,0,WACV)
                else:
                    flagSizer.Add(wx.StaticText(self.panel,label='na'),0,WACV)
            
        mainSizer.Add(flagSizer,0)
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        CancelBtn = wx.Button(self.panel,-1,'Cancel')
        CancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add(CancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
        
    def GetSelection(self):
        return self.newflags

    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_CANCEL)

################################################################################
def G2MessageBox(parent,msg,title='Error'):
    '''Simple code to display a error or warning message

    TODO: replace wx.MessageDialog with one derived from wx.Dialog because
    on most platforms wx.MessageDialog is a native widget and CentreOnParent
    will not function. 
    '''
    dlg = wx.MessageDialog(parent,StripIndents(msg), title, wx.OK|wx.CENTRE)
    dlg.CentreOnParent()
    dlg.ShowModal()
    dlg.Destroy()

def ShowScrolledInfo(parent,txt,width=600,height=400,header='Warning info',
                         buttonlist=None):
    '''Simple code to display possibly extensive error or warning text
    in a scrolled window.

    :param wx.Frame parent: parent window for 
    :param str txt: text to be displayed
    :param int width: lateral of window in pixels (defaults to 600)
    :param int height: vertical dimension of window in pixels (defaults to 400)
    :param str header: title to be placed on window
    :param list buttonlist: list of button Ids to show, or one or more 
      pairs of values, where the first is a label to place on the button 
      and the second is a routine that is called if the button is pressed.
      The default is None which places a single "Close" button that 
      returns wx.ID_CANCEL
    :returns: the wx Id for the selected button

    example::

       res = ShowScrolledInfo(self.frame,msg,header='Please Note',buttonlist=[
               ('Open', lambda event: event.GetEventObject().GetParent().EndModal(wx.ID_OK)),
               ('Skip', lambda event: event.GetEventObject().GetParent().EndModal(wx.ID_CANCEL))
              ])
       if res == wx.ID_OK:
           pass
    '''
    
    dlg = wx.Dialog(parent.GetTopLevelParent(),wx.ID_ANY,header, style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    spanel = wxscroll.ScrolledPanel(dlg, wx.ID_ANY, size=(width-20, height))
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(spanel,1,wx.ALL|wx.EXPAND,1)

    txtSizer = wx.BoxSizer(wx.VERTICAL)
    txt = wx.StaticText(spanel,wx.ID_ANY,txt)
    txt.Wrap(width-20)
    txt.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
    txtSizer.Add(txt,1,wx.ALL|wx.EXPAND,1)
    spanel.SetSizer(txtSizer)
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    if buttonlist is None:
        btn = wx.Button(dlg, wx.ID_CLOSE) 
        btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
        btnsizer.Add(btn)
    else:
        for b in buttonlist:
            if isinstance(b, (list, tuple)):
                btn = wx.Button(dlg, wx.ID_ANY, b[0]) 
                btn.Bind(wx.EVT_BUTTON,b[1])
            else:
                btn = wx.Button(dlg, b) 
                btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(event.Id))
            btnsizer.Add(btn)
    mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
    dlg.SetSizer(mainSizer)
    mainSizer.Fit(dlg)
    spanel.SetAutoLayout(1)
    spanel.SetupScrolling()
    #dlg.SetMaxSize((-1,400)) 
    dlg.CenterOnParent()
    ans = dlg.ShowModal()
    dlg.Destroy()
    return ans

def ShowScrolledColText(parent,txt,width=600,height=400,header='Warning info',col1len=999):
    '''Simple code to display tabular information in a scrolled wx.Dialog 
    window.

    Lines ending with a colon (:) are centered across all columns 
    and have a grey background. 
    Lines beginning and ending with '**' are also are centered 
    across all columns and are given a yellow background
    All other lines have columns split by tab (\\t) characters. 

    :param wx.Frame parent: parent window 
    :param str txt: text to be displayed
    :param int width: lateral of window in pixels (defaults to 600)
    :param int height: vertical dimension of window in pixels (defaults to 400)
    :param str header: title to be placed on window
    '''
    
    dlg = wx.Dialog(parent.GetTopLevelParent(),wx.ID_ANY,header, style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    spanel = wxscroll.ScrolledPanel(dlg, wx.ID_ANY, size=(width-20, height))
    spanel.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(spanel,1,wx.ALL|wx.EXPAND,1)

    cols = 1
    for i,line in enumerate(txt.split('\n')):
        cols = max(cols,line.count('\t')+1)

    txtSizer = wx.GridBagSizer(0,9)
    for i,line in enumerate(txt.split('\n')):
        if line.strip().endswith(':'):
            st = wx.StaticText(spanel,wx.ID_ANY,line)
            txtSizer.Add(st,pos=(i,0),span=(0,cols),flag=wx.EXPAND)
            st.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE))
            continue
        elif line.strip().startswith('**') and line.strip().endswith('**'):
            st = wx.StaticText(spanel,wx.ID_ANY,line,style=wx.ALIGN_CENTER)
            st.SetBackgroundColour(DULL_YELLOW)
            txtSizer.Add(st,pos=(i,0),span=(0,cols),flag=wx.EXPAND)
            continue
        items = line.split('\t')
        for col in range(cols):
            if col < len(items):
                item = items[col].strip()
            else:
                item = ''
            t = item[:]
            s = ''
            #if len(t) > col1len: GSASIIpath.IPyBreak()
            while col == 0 and len(t) > col1len:
                b = -1
                for sym in (') ',' * ',' + ',' - ',' && '):
                    b = max(b, t.rfind(sym,0,col1len))
                if b > 20:
                    s += t[:b+1] 
                    t = '\n\t' + t[b+1:]
                    continue
                break
            s += t
            st = wx.StaticText(spanel,wx.ID_ANY,s)
            if col == 0: st.Wrap(650)  # last resort...
            st.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            txtSizer.Add(st,pos=(i,col),flag=wx.EXPAND)
        #txtSizer.AddGrowableRow(i)
    txtSizer.AddGrowableCol(0)  #to fill screen
    spanel.SetSizer(txtSizer)
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = wx.Button(dlg, wx.ID_CLOSE)
    btn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_CANCEL))
    btnsizer.Add(btn)
    mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
    dlg.SetSizer(mainSizer)
    mainSizer.Fit(dlg)
    spanel.SetAutoLayout(1)
    spanel.SetupScrolling()
    #dlg.SetMaxSize((-1,400)) 
    dlg.CenterOnParent()
    dlg.ShowModal()
    dlg.Destroy()
    
def G2ScrolledGrid(G2frame,lbl,title,tbl,colLbls,colTypes,maxSize=(600,300),comment=''):
    '''Display a scrolled table of information in a dialog window

    :param wx.Frame G2frame: parent for dialog
    :param str lbl: label for window 
    :param str title: window title
    :param list tbl: list of lists where inner list is each row
    :param list colLbls: list of str with labels for each column
    :param list colTypes: Data types for each column (such as 
      wg.GRID_VALUE_STRING,wg.GRID_VALUE_FLOAT)
    :param list maxSize: Maximum size for the table in points. Defaults to 
      (600,300)
      :param str comment: optional text that appears below table

    Example::

       row = ['item1',1.234,'description of item']
       colTypes = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_FLOAT+':8,4',wg.GRID_VALUE_STRING]
       colLbls = ['item name','value','Description']
       G2ScrolledGrid(frm,'window label','title',20*[row],colLbls,colTypes)

    '''
        
    dlg = wx.Dialog(G2frame,title=title,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(wx.StaticText(dlg,label=lbl),
                      0,wx.ALIGN_CENTER_HORIZONTAL,0)
    sizer.Add((-1,15))
    rowlbl = [str(i+1) for i in range(len(tbl))]
    wxtbl = Table(tbl,rowLabels=rowlbl,colLabels=colLbls,types=colTypes)
        
    scrGrid = wx.ScrolledWindow(dlg)
    wxGrid = GSGrid(scrGrid)
    wxGrid.SetTable(wxtbl, True)
    wxGrid.AutoSizeColumns(False)
    wxGrid.EnableEditing(False)
    gridSizer = wx.BoxSizer(wx.VERTICAL)
    gridSizer.Add(wxGrid,1,wx.EXPAND,1)
    gridSizer.Layout()
    Size = gridSizer.GetMinSize()
    
    Size[0] = min(Size[0]+25,maxSize[0])
    Size[1] = min(Size[1]+25,maxSize[1])
    scrGrid.SetSizer(gridSizer)
    scrGrid.SetMinSize(Size)
    scrGrid.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1))
    scrGrid.Scroll(0,0)
    sizer.Add(scrGrid,1,wx.EXPAND,1)
    if len(comment):
        sizer.Add(wx.StaticText(dlg,label=comment))
        
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    btnsizer.Add((-1,-1),1,wx.EXPAND,1)
    btn = wx.Button(dlg, wx.ID_OK)
    btn.SetDefault()
    btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_OK))
    btnsizer.Add(btn)
    btnsizer.Add((-1,-1),1,wx.EXPAND,1)
    sizer.Add(btnsizer, 0, wx.EXPAND|wx.ALL, 5)
       
    sizer.Layout()    
    dlg.SetSizer(sizer)
    sizer.Fit(dlg)
    dlg.CenterOnParent()
    dlg.ShowModal()
    dlg.Destroy()

################################################################################
class PickTwoDialog(wx.Dialog):
    '''This does not seem to be in use
    '''
    def __init__(self,parent,title,prompt,names,choices):
        wx.Dialog.__init__(self,parent,-1,title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = None
        self.prompt = prompt
        self.choices = choices
        self.names = names
        self.Draw()

    def Draw(self):
        Indx = {}
        
        def OnSelection(event):
            Obj = event.GetEventObject()
            id = Indx[Obj.GetId()]
            self.choices[id] = Obj.GetValue().encode()  #to avoid Unicode versions
            self.Draw()
            
        if self.panel:
            self.panel.DestroyChildren()  #safe: wx.Panel
            self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,self.prompt),0,wx.ALIGN_CENTER)
        for isel,name in enumerate(self.choices):
            lineSizer = wx.BoxSizer(wx.HORIZONTAL)
            lineSizer.Add(wx.StaticText(self.panel,-1,'Reference atom '+str(isel+1)),0,wx.ALIGN_CENTER)
            nameList = self.names[:]
            if isel:
                if self.choices[0] in nameList:
                    nameList.remove(self.choices[0])
            choice = wx.ComboBox(self.panel,-1,value=name,choices=nameList,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            Indx[choice.GetId()] = isel
            choice.Bind(wx.EVT_COMBOBOX, OnSelection)
            lineSizer.Add(choice,0,WACV)
            mainSizer.Add(lineSizer)
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        CancelBtn = wx.Button(self.panel,-1,'Cancel')
        CancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add(CancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
        
    def GetSelection(self):
        return self.choices

    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_CANCEL)

################################################################################
class SingleFloatDialog(wx.Dialog):
    '''Dialog to obtain a single float value from user

    :param wx.Frame parent: name of parent frame
    :param str title: title string for dialog
    :param str prompt: string to tell user what they are inputing
    :param str value: default input value, if any
    :param list limits: upper and lower value used to set bounds for entry, use [None,None]
      for no bounds checking, [None,val] for only upper bounds, etc. Default is [0,1].
      Values outside of limits will be ignored.
    :param str format: string to format numbers. Defaults to '%.5g'. Use '%d' to have
      integer input (but dlg.GetValue will still return a float). 
    
    Typical usage::

            limits = (0,1)
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new value for...',default,limits)
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
            dlg.Destroy()    

    '''
    def __init__(self,parent,title,prompt,value,limits=[0.,1.],fmt='%.5g'):
        wx.Dialog.__init__(self,parent,-1,title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.CenterOnParent()
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self,-1,prompt),0,wx.ALIGN_CENTER)
        #valItem = wx.TextCtrl(self,-1,value=self.format%(self.value),style=wx.TE_PROCESS_ENTER)
        self.buffer = {0:float(fmt%(value))}
        a,b = fmt[1:].split('.')
        f = b[-1]
        try:
            d = int(b[:-1])
        except:
            d = 5
        try:
            w = int(a)
        except:
            w = 3+d
        self.OKbtn = wx.Button(self,wx.ID_OK)
        CancelBtn = wx.Button(self,wx.ID_CANCEL)
        valItem = ValidatedTxtCtrl(self,self.buffer,0,nDig=(w,d,f),
                                       xmin=limits[0],xmax=limits[1],
                                       OKcontrol=self.ControlOKButton)
        mainSizer.Add(valItem,0,wx.ALIGN_CENTER)
        self.OKbtn.Bind(wx.EVT_BUTTON, self.OnOk)
        CancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(self.OKbtn)
        btnSizer.Add(CancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.SetSizer(mainSizer)
        self.Fit()

    def GetValue(self):
        return self.buffer[0]
        
    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_CANCEL)

    def ControlOKButton(self,setvalue):
        '''Enable or Disable the OK button for the dialog. Note that this is
        passed into the ValidatedTxtCtrl for use by validators.

        :param bool setvalue: if True, all entries in the dialog are
          checked for validity. if False then the OK button is disabled.

        '''
        if setvalue: # turn button on, do only if all controls show as valid
            self.OKbtn.Enable()
        else:
            self.OKbtn.Disable()

class SingleIntDialog(SingleFloatDialog):
    '''Dialog to obtain a single int value from user

    :param wx.Frame parent: name of parent frame
    :param str title: title string for dialog
    :param str prompt: string to tell user what they are inputing
    :param str value: default input value, if any
    :param list limits: upper and lower value used to set bounds for entries. Default
      is [None,None] -- for no bounds checking; use [None,val] for only upper bounds, etc.
      Default is [0,1]. Values outside of limits will be ignored.
    
    Typical usage::

            limits = (0,None)  # allows zero or positive values only
            dlg = G2G.SingleIntDialog(G2frame,'New value','Enter new value for...',default,limits)
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
            dlg.Destroy()    

    '''
    def __init__(self,parent,title,prompt,value,limits=[None,None]):
        SingleFloatDialog.__init__(self,parent,title,prompt,value,limits=limits,format='%d')
    def GetValue(self):
        return int(self.value)

################################################################################
class MultiDataDialog(wx.Dialog):
    '''Dialog to obtain multiple values from user. Use ``dlg.GetValues()`` to 
    get the values set in the window.

    :param wx.Frame parent: parent frame for dialog to be created
    :param str title: title to place on top of window
    :param list prompts: a string to describe each item. Each entry
      in this list will designate a row in the generated window. 
    :param list values: a list of initial values for each item. Use a 
      nested list when multiple entries are placed on a single row of 
      the window (see discussion of formats, below). Number of items
      in the outer list should match the length of ``prompts``. The 
      total number of items should match ``formats``.
    :param list limits: A nested list with an upper and lower value 
       for each item or for a choice/edit control a list of allowed 
       values. Use a nested list when multiple entries are placed on 
       a single row of the window (see discussion of formats, below).
       Number of items in the outer list should match the length of 
       ``prompts``. The total number of items should match ``formats``.
    :param list testfxns: A nested list of string test functions.
       The total number of items should match ``formats`` or should be
       left as the default (None).
    :param list formats: A list of values for each entry in the 
       window. Several different types of values are possible: 

       * An "old-style" format string (e.g. ``%5d`` or ``%.3f``) 
         which will be used to display each item's value

       * Or a keyword that specifies how the values are used. 
         Allowed keywords are:

         * ``choice``: for a pull-down list;
         * ``bool``: for a yes/no checkbox;
         * ``str``: for a text entry 
         * ``edit``: for a pull-down list that allows one to enter an arbitrary value.

       * Alternately, a value can be a list of items, in which case multiple 
         entries are placed on a single row of the window. When this is done, 
         any value in the list other than ``choice`` or ``edit`` is used as 
         text to be placed between the ComboBoxes.

       The number of items in the outer list should match the length 
       of ``prompts``.
    :param str header: a string to be placed at the top of the 
       window. Ignored if None (the default.)

    example 1::

        dlg = G2G.MultiDataDialog(G2frame,title='ISOCIF search',
                prompts=['lattice constants tolerance',
                         'coordinate tolerance',
                         'occupancy tolerance'],
                values=[0.001,0.01,0.1],
                limits=3*[[0.,2.]],formats=3*['%.5g'],
                header=isoCite)
        dlg.ShowModal()
        latTol,coordTol,occTol = dlg.GetValues()
        dlg.Destroy()

    example 2::

        nm = [' ','0','1','-1','2','-2','3','-3','4','5','6','7','8','9']
        dm = ['1','2','3','4','5','6']
        kfmt = ['choice','/','choice',',    ','choice','/','choice',',    ','choice','/','choice',' ']
        dlg = MultiDataDialog(G2frame,title='options',
                prompts=[' k-vector 1 (x,y,z)',
                         ' k-vector 2 (x,y,z)'],
                values=[3*['0','','2',''],3*[' ','','2','']],
                limits=[3*[nm[1:],'',dm,''],3*[nm,'',dm,'']],
                formats=[kfmt,kfmt])
        if dlg.ShowModal() == wx.ID_OK: print(dlg.GetValues())
'''
    def __init__(self,parent,title,prompts,values,limits=[[0.,1.],],
                     testfxns=None,formats=['%.5g',],header=None):
        wx.Dialog.__init__(self,parent,-1,title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = None
        self.limits = limits
        self.values = values
        self.prompts = prompts
        self.formats = formats
        if testfxns is None:
            self.testfxns = []
            for i in formats:
                if type(i) is list:
                    self.testfxns.append(len(i)*[None])
                else:
                    self.testfxns.append(None)                    
        else:
            self.testfxns = testfxns
        self.header = header
        self.Draw()
        
    def Draw(self):
        
        def OnEditItem(event):
            if event: event.Skip()
            Obj = event.GetEventObject()
            fmt = Indx[Obj][-1]
            if type(fmt) is list:
                tid,idl,limits = Indx[Obj][:3]
            else:
                tid,idl = Indx[Obj],0
            val = Obj.GetValue()
            try:
                eval(val)
                self.values[tid][idl] = val
                return
            except:
                pass
            try:   # deal with mixed fractions (example: 1 3/4)
                val = val.replace(' ','+')
                val = str(eval(val))
                self.values[tid][idl] = val
                return
            except:
                pass
                    
        def OnValItem(event):
            if event: event.Skip()
            Obj = event.GetEventObject()
            fmt = Indx[Obj][-1]
            if type(fmt) is list:
                tid,idl,limits = Indx[Obj][:3]
                if 'testfxn' in fmt:
                    testfxn = Indx[Obj][3][tid]
                    val = Obj.GetValue().strip()
                    if testfxn(val):
                        self.values[tid][idl] = val
                        Obj.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE))
                        Obj.SetForegroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT))
                    else:
                        Obj.SetBackgroundColour(wx.YELLOW)
                        Obj.SetForegroundColour("red")
                    Obj.SetValue('%s'%(val))                    
                self.values[tid][idl] = Obj.GetValue()
            elif 'bool' in fmt:
                self.values[Indx[Obj][0]] = Obj.GetValue()
            elif 'str' in fmt:
                tid,limits = Indx[Obj][:2]
                try:
                    val = Obj.GetValue()
                    if val not in limits:
                        raise ValueError
                except ValueError:
                    val = self.values[tid]
                self.values[tid] = val
                Obj.SetValue('%s'%(val))
            elif 'testfxn' in fmt:
                tid,x,testfxn = Indx[Obj][:3]
                val = Obj.GetValue()
                if testfxn(val):
                    self.values[tid] = val
                    Obj.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE))
                    Obj.SetForegroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNTEXT))
                else:
                    Obj.SetBackgroundColour(wx.YELLOW)
                    Obj.SetForegroundColour("red")
                Obj.SetValue('%s'%(val))                    
            elif 'choice' in fmt:
                self.values[Indx[Obj][0]] = Obj.GetValue()
            
        Indx = {}
        if self.panel: self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if self.header:
            txt = wx.StaticText(self.panel,wx.ID_ANY,self.header)
            txt.Wrap(400)
            mainSizer.Add(txt)
            mainSizer.Add((-1,5))
            HorizontalLine(mainSizer,self.panel)
            mainSizer.Add((-1,5))
        lineSizer = wx.FlexGridSizer(0,2,5,5)
        for tid,[prompt,value,limits,testfxn,fmt] in enumerate(zip(self.prompts,self.values,self.limits,self.testfxns,self.formats)):
            lineSizer.Add(wx.StaticText(self.panel,label=prompt),0,wx.ALIGN_CENTER)
            if type(fmt) is list:
                valItem = wx.BoxSizer(wx.HORIZONTAL)
                for idl,item in enumerate(fmt):
                    if value[idl] in limits[idl]:
                        initVal = value[idl]
                    else:
                        initVal = limits[idl][0]
                    if 'edit' in item:
                        style = wx.CB_DROPDOWN
                        listItem = wx.ComboBox(self.panel,value=initVal,
                            choices=limits[idl],style=style)
                        Indx[listItem] = [tid,idl,limits,testfxn,fmt]
                        listItem.Bind(wx.EVT_TEXT,OnEditItem)
                        listItem.Bind(wx.EVT_COMBOBOX,OnValItem)
                    elif 'testfxn' in item:
                        listItem = wx.TextCtrl(self.panel,value='%s'%(value),style=wx.TE_PROCESS_ENTER)
                        Indx[listItem] = [tid,idl,limits,testfxn,fmt]
                        listItem.Bind(wx.EVT_TEXT_ENTER,OnValItem)
                        listItem.Bind(wx.EVT_KILL_FOCUS,OnValItem)
                        listItem.SetValue('%s'%value[idl])
                    elif 'choice' in item:
                        style = wx.CB_READONLY|wx.CB_DROPDOWN
                        listItem = wx.ComboBox(self.panel,value=initVal,
                            choices=limits[idl],style=style)
                        Indx[listItem] = [tid,idl,limits,testfxn,fmt]
                        listItem.Bind(wx.EVT_COMBOBOX,OnValItem)
                    else:
                        listItem = wx.StaticText(self.panel,wx.ID_ANY,item)
                    valItem.Add(listItem,0,WACV)
            elif 'bool' in fmt:
                valItem = wx.CheckBox(self.panel,label='')
                valItem.Bind(wx.EVT_CHECKBOX,OnValItem)
                valItem.SetValue(value)
            elif 'str' in fmt:
                valItem = wx.TextCtrl(self.panel,value='%s'%(value),style=wx.TE_PROCESS_ENTER)
                valItem.Bind(wx.EVT_TEXT_ENTER,OnValItem)
                valItem.Bind(wx.EVT_KILL_FOCUS,OnValItem)
                valItem.SetValue('%s'%value)
            elif 'edit' in fmt:
                valItem = wx.ComboBox(self.panel,value=limits[0],choices=limits,style=wx.CB_DROPDOWN)
                valItem.Bind(wx.EVT_COMBOBOX,OnValItem)
                valItem.Bind(wx.EVT_TEXT,OnEditItem)
            elif 'choice' in fmt:
                valItem = wx.ComboBox(self.panel,value=limits[0],choices=limits,style=wx.CB_READONLY|wx.CB_DROPDOWN)
                valItem.Bind(wx.EVT_COMBOBOX,OnValItem)
            elif 'testfxn' in fmt:
                valItem = wx.TextCtrl(self.panel,value='%s'%(value),style=wx.TE_PROCESS_ENTER)
                valItem.Bind(wx.EVT_TEXT_ENTER,OnValItem)
                valItem.Bind(wx.EVT_KILL_FOCUS,OnValItem)
            else:
                if '%' in fmt:
                    if 'd' in fmt:
                        nDig = None
                    else:
                        sfmt = fmt[1:].split('.')
                        if not sfmt[0]: sfmt[0] = '10'
                        nDig = (int(sfmt[0]),int(sfmt[1][:-1]),sfmt[1][-1])
                valItem = ValidatedTxtCtrl(self.panel,self.values,tid,nDig=nDig,xmin=limits[0],xmax=limits[1])
                # valItem = wx.TextCtrl(self.panel,value=fmt%(value),style=wx.TE_PROCESS_ENTER)
                # valItem.Bind(wx.EVT_TEXT_ENTER,OnValItem)
                # valItem.Bind(wx.EVT_KILL_FOCUS,OnValItem)
            if type(fmt) is not list:
                Indx[valItem] = [tid,limits,fmt]
            lineSizer.Add(valItem,0,wx.ALIGN_CENTER)
        mainSizer.Add(lineSizer)
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        CancelBtn = wx.Button(self.panel,-1,'Cancel')
        CancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        OkBtn.SetDefault()
        btnSizer.Add(CancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()

    def GetValues(self):
        return self.values
        
    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_CANCEL)

################################################################################
class SingleStringDialog(wx.Dialog):
    '''Dialog to obtain a single string value from user
    
    :param wx.Frame parent: name of parent frame
    :param str title: title string for dialog
    :param str prompt: string to tell use what they are inputting
    :param str value: default input value, if any
    :param tuple size: specifies default size and width for the text 
      entry section of the dialog [default (200,-1)]. If the vertical 
      size (the second number) is greater than 20 (~ a single line) then 
      the textbox will allow inclusion of new-line characters. In single-line
      mode, return causes the dialog to close.
    :param str help: if supplied, a help button is added to the dialog that
      can be used to display the supplied help text/URL for setting this 
      variable. (Default is '', which is ignored.)
    :param list choices: a set of strings that provide optional values that 
      can be selected from; these can be edited if desired.
    '''
    def __init__(self,parent,title,prompt,value='',size=(200,-1),help='',
                     choices=None):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title,pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.value = value
        self.prompt = prompt
        self.CenterOnParent()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer1.Add(wx.StaticText(self.panel,-1,self.prompt),0,wx.ALIGN_CENTER)
        sizer1.Add((-1,-1),1,wx.EXPAND)
        if help:
            sizer1.Add(HelpButton(self.panel,help),0,wx.ALL)
        mainSizer.Add(sizer1,0,wx.EXPAND)
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        if choices:
            self.valItem = wx.ComboBox(self.panel, wx.ID_ANY, value,
                                size=size,choices=[value]+choices,
                                style=wx.CB_DROPDOWN|wx.TE_PROCESS_ENTER)
        elif size[1] > 20:
            self.valItem = wx.TextCtrl(self.panel,-1,value=self.value,size=size,
                                style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER)
        else:
            self.valItem = wx.TextCtrl(self.panel,-1,value=self.value,size=size)
        #sizer1.Add(self.valItem,0,wx.ALIGN_CENTER)
        #mainSizer.Add(sizer1,0,wx.EXPAND)
        mainSizer.Add(self.valItem,1,wx.EXPAND)
        btnsizer = wx.StdDialogButtonSizer()
        OKbtn = wx.Button(self.panel, wx.ID_OK)
        OKbtn.SetDefault()
        btnsizer.AddButton(OKbtn)
        btn = wx.Button(self.panel, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()

    def Show(self):
        '''Use this method after creating the dialog to post it
        :returns: True if the user pressed OK; False if the User pressed Cancel
        '''
        if self.ShowModal() == wx.ID_OK:
            self.value = self.valItem.GetValue()
            return True
        else:
            return False

    def GetValue(self):
        '''Use this method to get the value entered by the user
        :returns: string entered by user
        '''
        return self.value

################################################################################
class MultiStringDialog(wx.Dialog):
    '''Dialog to obtain a multi string values from user
    
    :param wx.Frame parent: name of parent frame
    :param str title: title string for dialog
    :param list prompts: list of strings to tell user what they are inputting
    :param list values: list of str default input values, if any
    :param int size: length of the input box in pixels
    :param bool addRows: if True, users can add rows to the table 
      (default is False)
    :param str hlp: if supplied, a help button is added to the dialog that
      can be used to display the supplied help text in this variable.
    :param str lbl: label placed at top of dialog  
    :returns: a wx.Dialog instance
    '''
    def __init__(self,parent,title,prompts,values=[],size=-1,
                     addRows=False,hlp=None, lbl=None):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title, 
                           pos=wx.DefaultPosition,
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.values = list(values)
        self.prompts = list(prompts)
        self.addRows = addRows
        self.size = size
        self.hlp = hlp
        self.lbl = lbl
        self.CenterOnParent()
        self.Paint()

    def Paint(self):
        if self.GetSizer():
            self.GetSizer().Clear(True)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if self.hlp:
            btnsizer = wx.BoxSizer(wx.HORIZONTAL)
            hlp = HelpButton(self, self.hlp, wrap=450)
            btnsizer.Add((-1,-1),1, wx.EXPAND, 1)
            #btnsizer.Add(hlp,0,wx.ALIGN_RIGHT|wx.ALL)
            btnsizer.Add(hlp,0)
            mainSizer.Add(btnsizer,0,wx.EXPAND)
        if self.lbl:
            mainSizer.Add(wx.StaticText(self,wx.ID_ANY,self.lbl))
            mainSizer.Add((-1,15))
        promptSizer = wx.FlexGridSizer(0,2,5,5)
        promptSizer.AddGrowableCol(1,1)
        self.Indx = {}
        for prompt,value in zip(self.prompts,self.values):
            promptSizer.Add(wx.StaticText(self,-1,prompt))
            valItem = wx.TextCtrl(self,-1,value=value,style=wx.TE_PROCESS_ENTER,size=(self.size,-1))
            self.Indx[valItem.GetId()] = prompt
            valItem.Bind(wx.EVT_TEXT,self.newValue)
            promptSizer.Add(valItem,1,wx.EXPAND,1)
        mainSizer.Add(promptSizer,1,wx.ALL|wx.EXPAND,1)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        OKbtn = wx.Button(self, wx.ID_OK)
        OKbtn.SetDefault()
        btnsizer.Add((1,1),1,wx.EXPAND,1)
        btnsizer.Add(OKbtn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.Add(btn)
        btnsizer.Add((1,1),1,wx.EXPAND,1)
        if self.addRows:
            btn = wx.Button(self, wx.ID_ANY,'+',style=wx.BU_EXACTFIT)
            btn.Bind(wx.EVT_BUTTON,self.onExpand)
            btnsizer.Add(btn)
        mainSizer.Add(btnsizer,0,wx.EXPAND)
        self.SetSizer(mainSizer)
        self.Fit()

    def onExpand(self,event):
        self.values.append('')
        self.prompts.append('item '+str(len(self.values)))
        self.Paint()
        
    def newValue(self,event):
        Obj = event.GetEventObject()
        item = self.Indx[Obj.GetId()]
        id = self.prompts.index(item)
        self.values[id] = Obj.GetValue()

    def Show(self):
        '''Use this method after creating the dialog to post it
        
        :returns: True if the user pressed OK; False if the User pressed Cancel
        '''
        if self.ShowModal() == wx.ID_OK:
            return True
        else:
            return False

    def GetValues(self):
        '''Use this method to get the value(s) entered by the user
        
        :returns: a list of strings entered by user
        '''
        return self.values

################################################################################
class G2ColumnIDDialog(wx.Dialog):
    '''A dialog for matching column data to desired items; some columns may be ignored.
    
    :param wx.Frame ParentFrame: reference to parent frame
    :param str title: heading above list of choices
    :param str header: Title to place on window frame 
    :param list ChoiceList: a list of possible choices for the columns
    :param list ColumnData: lists of column data to be matched with ChoiceList
    :param bool monoFont: If False (default), use a variable-spaced font;
      if True use a equally-spaced font.
    :param kw: optional keyword parameters for the wx.Dialog may
      be included such as size [which defaults to `(320,310)`] and
      style (which defaults to ``wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.CENTRE | wx.OK | wx.CANCEL``);
      note that ``wx.OK`` and ``wx.CANCEL`` controls
      the presence of the eponymous buttons in the dialog.
    :returns: the name of the created dialog
    
    '''

    def __init__(self,parent, title, header,Comments,ChoiceList, ColumnData,
                 monoFont=False, **kw):

        def OnOk(sevent):
            OK = True
            selCols = []
            for col in self.sel:
                item = col.GetValue()
                if item != ' ' and item in selCols:
                    OK = False
                    break
                else:
                    selCols.append(item)
            parent = self.GetParent()
            if not OK:
                parent.ErrorDialog('Duplicate',item+' selected more than once')
                return
            if parent is not None: parent.Raise()
            self.EndModal(wx.ID_OK)
            
        def OnModify(event):
            if event: event.Skip()
            Obj = event.GetEventObject()
            icol,colData = Indx[Obj.GetId()]
            modify = Obj.GetValue()
            if not modify:
                return
            #print 'Modify column',icol,' by', modify
            for i,item in enumerate(self.ColumnData[icol]):
                self.ColumnData[icol][i] = str(eval(item+modify))
            colData.SetValue('\n'.join(self.ColumnData[icol]))
            Obj.SetValue('')
            
        # process keyword parameters, notably style
        options = {'size':(600,310), # default Frame keywords
                   'style':wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.CENTRE| wx.OK | wx.CANCEL,
                   }
        options.update(kw)
        self.Comments = ''.join(Comments)
        self.ChoiceList = ChoiceList
        self.ColumnData = ColumnData
        if options['style'] & wx.OK:
            useOK = True
            options['style'] ^= wx.OK
        else:
            useOK = False
        if options['style'] & wx.CANCEL:
            useCANCEL = True
            options['style'] ^= wx.CANCEL
        else:
            useCANCEL = False        
        # create the dialog frame
        wx.Dialog.__init__(self,parent,wx.ID_ANY,header,**options)
        panel = wxscroll.ScrolledPanel(self)
        # fill the dialog
        Sizer = wx.BoxSizer(wx.VERTICAL)
        Sizer.Add((-1,5))
        if self.Comments:
            if title[-1] == ':':
                title = title[:-1] + ' using header line(s):'
            else:
                title += ' using header line(s):'
            Sizer.Add(wx.StaticText(panel,label=title),0)
            Sizer.Add((5,5))
            if self.Comments[-1] != '\n': self.Comments += '\n'
            txt = wx.StaticText(panel,label=self.Comments)
            txt.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            font1 = wx.Font(txt.GetFont().GetPointSize(),wx.MODERN, wx.NORMAL, wx.NORMAL, False)
            txt.SetFont(font1)
            Sizer.Add(txt,0,wx.ALL|wx.EXPAND,0)
            txtSize = txt.GetSize()[1]
        else:
            Sizer.Add(wx.StaticText(panel,label=title),0)
            txtSize = 0
        columnsSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sel = []
        self.mod = []
        Indx = {}
        for icol,col in enumerate(self.ColumnData):
            colSizer = wx.BoxSizer(wx.VERTICAL)
            colSizer.Add(wx.StaticText(panel,label=' Column #%d Select:'%(icol)))
            self.sel.append(wx.ComboBox(panel,value=' ',choices=self.ChoiceList,style=wx.CB_READONLY|wx.CB_DROPDOWN))
            colSizer.Add(self.sel[-1])
            colData = wx.TextCtrl(panel,value='\n'.join(self.ColumnData[icol]),size=(120,-1),
                style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_DONTWRAP)
            colSizer.Add(colData,1,wx.ALL|wx.EXPAND,1)
            colSizer.Add(wx.StaticText(panel,label=' Modify by:'))
            mod = wx.TextCtrl(panel,size=(120,-1),value='',style=wx.TE_PROCESS_ENTER)
            mod.Bind(wx.EVT_TEXT_ENTER,OnModify)
            mod.Bind(wx.EVT_KILL_FOCUS,OnModify)
            Indx[mod.GetId()] = [icol,colData]
            colSizer.Add(mod)
            columnsSizer.Add(colSizer,0,wx.ALL|wx.EXPAND,10)
        Sizer.Add(columnsSizer,1,wx.ALL|wx.EXPAND,1)
        Sizer.Add(wx.StaticText(panel,label=' For modify by, enter arithmetic string eg. "-12345.67". "+", "-", "*", "/", "**" all allowed'),0) 
        Sizer.Add((-1,10))
        # OK/Cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        if useOK:
            self.OKbtn = wx.Button(panel, wx.ID_OK)
            self.OKbtn.SetDefault()
            btnsizer.AddButton(self.OKbtn)
            self.OKbtn.Bind(wx.EVT_BUTTON, OnOk)
        if useCANCEL:
            btn = wx.Button(panel, wx.ID_CANCEL)
            btnsizer.AddButton(btn)
        btnsizer.Realize()
        Sizer.Add((-1,5))
        Sizer.Add(btnsizer,0,wx.ALIGN_LEFT,20)
        Sizer.Add((-1,5))
        # OK done, let's get outa here
        panel.SetSizer(Sizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
        Size = [450,375]
        panel.SetSize(Size)
        Size[0] += 25; Size[1]+= 25+txtSize
        self.SetSize(Size)
        
    def GetSelection(self):
        'Returns the selected sample parm for each column'
        selCols = []
        for item in self.sel:
            selCols.append(item.GetValue())
        return selCols,self.ColumnData
    
################################################################################
class G2HistoDataDialog(wx.Dialog):
    '''A dialog for editing histogram data globally.
    
    :param wx.Frame ParentFrame: reference to parent frame
    :param str title: heading above list of choices
    :param str header: Title to place on window frame 
    :param list ParmList: a list of names for the columns
    :param list ParmFmt: a list of formatting strings for the columns
    :param list: HistoList: a list of histogram names
    :param list ParmData: a list of lists of data matched to ParmList; one for each item in HistoList
    :param bool monoFont: If False (default), use a variable-spaced font;
      if True use a equally-spaced font.
    :param kw: optional keyword parameters for the wx.Dialog may
      be included such as size [which defaults to `(320,310)`] and
      style (which defaults to 
      ``wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.CENTRE | wx.OK | wx.CANCEL``);
      note that ``wx.OK`` and ``wx.CANCEL`` controls the presence of the eponymous buttons in the dialog.
    :returns: the modified ParmData
    
    '''

    def __init__(self,parent, title, header,ParmList,ParmFmt,HistoList,ParmData,
                 monoFont=False, **kw):

        def OnOk(sevent):
            if parent is not None: parent.Raise()
            self.EndModal(wx.ID_OK)
            
        def OnModify(event):
            Obj = event.GetEventObject()
            irow,it = Indx[Obj.GetId()]
            try:
                val = float(Obj.GetValue())
            except ValueError:
                val = self.ParmData[irow][it]
            self.ParmData[irow][it] = val
            Obj.SetValue(self.ParmFmt[it]%val)
                        
        # process keyword parameters, notably style
        options = {'size':(600,310), # default Frame keywords
                   'style':wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.CENTRE| wx.OK | wx.CANCEL,
                   }
        options.update(kw)
        self.ParmList = ParmList
        self.ParmFmt = ParmFmt
        self.HistoList = HistoList
        self.ParmData = ParmData
        nCol = len(ParmList)
        if options['style'] & wx.OK:
            useOK = True
            options['style'] ^= wx.OK
        else:
            useOK = False
        if options['style'] & wx.CANCEL:
            useCANCEL = True
            options['style'] ^= wx.CANCEL
        else:
            useCANCEL = False        
        # create the dialog frame
        wx.Dialog.__init__(self,parent,wx.ID_ANY,header,**options)
        panel = wxscroll.ScrolledPanel(self)
        # fill the dialog
        Sizer = wx.BoxSizer(wx.VERTICAL)
        Sizer.Add((-1,5))
        Sizer.Add(wx.StaticText(panel,label=title),0)
        dataSizer = wx.FlexGridSizer(0,nCol+1,0,0)
        self.sel = []
        self.mod = []
        Indx = {}
        for item in ['Histogram',]+self.ParmList:
            dataSizer.Add(wx.StaticText(panel,-1,label=' %10s '%(item)),0,WACV)
        for irow,name in enumerate(self.HistoList):
            dataSizer.Add(wx.StaticText(panel,label=name),0,WACV|wx.LEFT|wx.RIGHT,10)
            for it,item in enumerate(self.ParmData[irow]):
                dat = wx.TextCtrl(panel,-1,value=self.ParmFmt[it]%(item),style=wx.TE_PROCESS_ENTER)
                dataSizer.Add(dat,0,WACV)
                dat.Bind(wx.EVT_TEXT_ENTER,OnModify)
                dat.Bind(wx.EVT_KILL_FOCUS,OnModify)
                Indx[dat.GetId()] = [irow,it]
        Sizer.Add(dataSizer)
        Sizer.Add((-1,10))
        # OK/Cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        if useOK:
            self.OKbtn = wx.Button(panel, wx.ID_OK)
            self.OKbtn.SetDefault()
            btnsizer.AddButton(self.OKbtn)
            self.OKbtn.Bind(wx.EVT_BUTTON, OnOk)
        if useCANCEL:
            btn = wx.Button(panel, wx.ID_CANCEL)
            btnsizer.AddButton(btn)
        btnsizer.Realize()
        Sizer.Add((-1,5))
        Sizer.Add(btnsizer,0,wx.ALIGN_LEFT,20)
        Sizer.Add((-1,5))
        # OK done, let's get outa here
        panel.SetSizer(Sizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
        Size = [450,375]
        panel.SetSize(Size)
        Size[0] += 25; Size[1]+= 25
        self.SetSize(Size)
        
    def GetData(self):
        'Returns the modified ParmData'
        return self.ParmData
    
################################################################################
def ItemSelector(ChoiceList, ParentFrame=None,
                 title='Select an item',
                 size=None, header='Item Selector',
                 useCancel=True,multiple=False):
    ''' Provide a wx dialog to select a single item or multiple items from list of choices

    :param list ChoiceList: a list of choices where one will be selected
    :param wx.Frame ParentFrame: Name of parent frame (default None)
    :param str title: heading above list of choices (default 'Select an item')
    :param wx.Size size: Size for dialog to be created (default None -- size as needed)
    :param str header: Title to place on window frame (default 'Item Selector')
    :param bool useCancel: If True (default) both the OK and Cancel buttons are offered
    :param bool multiple: If True then multiple items can be selected (default False)
    
    :returns: the selection index or None or a selection list if multiple is true

    Called by GSASIIdataGUI.OnReOrgSelSeq() Which is not fully implemented. 
    '''
    if multiple:
        if useCancel:
            dlg = G2MultiChoiceDialog(
                ParentFrame,title, header, ChoiceList)
        else:
            dlg = G2MultiChoiceDialog(
                ParentFrame,title, header, ChoiceList,
                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.OK|wx.CENTRE)
    else:
        if useCancel:
            dlg = wx.SingleChoiceDialog(
                ParentFrame,title, header, ChoiceList)
        else:
            dlg = wx.SingleChoiceDialog(
                ParentFrame,title, header,ChoiceList,
                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.OK|wx.CENTRE)
    if size: dlg.SetSize(size)
    if dlg.ShowModal() == wx.ID_OK:
        if multiple:
            dlg.Destroy()
            return dlg.GetSelections()
        else:
            dlg.Destroy()
            return dlg.GetSelection()
    else:
        dlg.Destroy()
        return None
    dlg.Destroy()

########################################################
# Column-order selection dialog
def GetItemOrder(parent,keylist,vallookup,posdict):
    '''Creates a dialog where items can be ordered into columns
    
    :param list keylist: is a list of keys for column assignments
    :param dict vallookup: is a dict keyed by names in keylist where each item is a dict. 
       Each inner dict contains variable names as keys and their associated values 
    :param dict posdict: is a dict keyed by names in keylist where each item is a dict. 
       Each inner dict contains column numbers as keys and their associated
       variable name as a value. This is used for both input and output.
       
    '''
    dlg = wx.Dialog(parent,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    sizer = wx.BoxSizer(wx.VERTICAL)
    spanel = OrderBox(dlg,keylist,vallookup,posdict)
    spanel.Fit()
    sizer.Add(spanel,1,wx.EXPAND)
    btnsizer = wx.StdDialogButtonSizer()
    btn = wx.Button(dlg, wx.ID_OK)
    btn.SetDefault()
    btnsizer.AddButton(btn)
    #btn = wx.Button(dlg, wx.ID_CANCEL)
    #btnsizer.AddButton(btn)
    btnsizer.Realize()
    sizer.Add(btnsizer, 0, wx.EXPAND|wx.ALL, 5)
    dlg.SetSizer(sizer)
    sizer.Fit(dlg)
    dlg.ShowModal()

################################################################################
class MultiIntegerDialog(wx.Dialog):
    '''Input a series of integers based on prompts
    '''
    def __init__(self,parent,title,prompts,values):
        wx.Dialog.__init__(self,parent,-1,title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.values = values
        self.prompts = prompts
        self.Draw()
        
    def Draw(self):
        
        def OnValItem(event):
            event.Skip()
            Obj = event.GetEventObject()
            ind = Indx[Obj.GetId()]
            try:
                val = int(Obj.GetValue())
                if val <= 0:
                    raise ValueError
            except ValueError:
                val = self.values[ind]
            self.values[ind] = val
            Obj.SetValue('%d'%(val))
            
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        Indx = {}
        for ival,[prompt,value] in enumerate(zip(self.prompts,self.values)):
            mainSizer.Add(wx.StaticText(self.panel,-1,prompt),0,wx.ALIGN_CENTER)
            valItem = wx.TextCtrl(self.panel,-1,value='%d'%(value),style=wx.TE_PROCESS_ENTER)
            mainSizer.Add(valItem,0,wx.ALIGN_CENTER)
            Indx[valItem.GetId()] = ival
            valItem.Bind(wx.EVT_TEXT_ENTER,OnValItem)
            valItem.Bind(wx.EVT_KILL_FOCUS,OnValItem)
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        CancelBtn = wx.Button(self.panel,-1,'Cancel')
        CancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add(CancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()

    def GetValues(self):
        return self.values
        
    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_CANCEL)

################################################################################
class MultiColumnSelection(wx.Dialog):
    '''Defines a Dialog widget that can be used to select an item from a multicolumn list.
    The first column should be short, but remaining columns are word-wrapped if the
    length of the information extends beyond the column. 
    
    When created, the dialog will be shown and <dlg>.Selection will be set to the index
    of the selected row, or -1. Be sure to use <dlg>.Destroy() to remove the window
    after reading the selection. If the dialog cannot be shown because a very old
    version of wxPython is in use, <dlg>.Selection will be None.

    If checkLbl is provided with a value, then a set of check buttons starts the table
    and <dlg>.Selections has the checked rows.
    
    :param wx.Frame parent: the parent frame (or None)
    :param str title: A title for the dialog window
    :param list colLabels: labels for each column
    :param list choices: a nested list with a value for each row in the table. Within each value
      should be a list of values for each column. There must be at least one value, but it is
      OK to have more or fewer values than there are column labels (colLabels). Extra are ignored
      and unspecified columns are left blank.
    :param list colWidths: a list of int values specifying the column width for each
      column in the table (pixels). There must be a value for every column label (colLabels).
    :param str checkLbl: A label for a row of checkboxes added at the beginning of the table.
       This option seems to be broken.
    :param int height: an optional height (pixels) for the table (defaults to 400)
    :param bool centerCols: if True, items in each column are centered. Default is False
    
    Example use::
    
        lbls = ('col 1','col 2','col 3')
        choices=(['test1','explanation of test 1'],
                 ['b', 'a really really long line that will be word-wrapped'],
                 ['test3','more explanation text','optional 3rd column text'])
        colWidths=[200,400,100]
        dlg = MultiColumnSelection(frm,'select tutorial',lbls,choices,colWidths)
        value = choices[dlg.Selection][0]
        dlg.Destroy()
    
    '''
    def __init__(self, parent, title, colLabels, choices, colWidths, checkLbl="",
                     height=400, centerCols=False, *args, **kw):
        if len(colLabels) != len(colWidths):
            raise ValueError('Length of colLabels) != colWidths')
        sizex = 20 # extra room for borders, etc.
        for i in colWidths: sizex += i
        wx.Dialog.__init__(self, parent, wx.ID_ANY, title, *args,
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER,
                           size=(sizex,height), **kw)
        self.Selections = len(choices)*[False]
        try:
            from wx.lib.wordwrap import wordwrap
            import wx.lib.agw.ultimatelistctrl as ULC
        except ImportError:
            self.Selection = None
            return
        self.Selection = -1
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.list = ULC.UltimateListCtrl(self, agwStyle=ULC.ULC_REPORT|ULC.ULC_HAS_VARIABLE_ROW_HEIGHT
                                         |ULC.ULC_HRULES|ULC.ULC_HRULES|ULC.ULC_SINGLE_SEL)
        if centerCols:
            colPosition = ULC.ULC_FORMAT_CENTER
        else:
            colPosition = ULC.ULC_FORMAT_LEFT
            
        if checkLbl:
            self.list.InsertColumn(0, checkLbl, width=8*len(checkLbl), format=colPosition)
            inc = 1
        else:
            inc = 0
        for i,(lbl,wid) in enumerate(zip(colLabels, colWidths)):
            self.list.InsertColumn(i+inc, lbl, width=wid, format=colPosition)
        for i,item in enumerate(choices):
            if item[0].startswith('   '):
                item[0] = '--- '+item[0].strip()
            if checkLbl:
                def OnCheck(event,row=i):
                    self.Selections[row] = event.EventObject.GetValue()
                c = wx.CheckBox(self.list)
                c.Bind(wx.EVT_CHECKBOX,OnCheck)
                self.list.InsertStringItem(i, "")
                citem = self.list.GetItem(i,0)
                citem.SetWindow(c)
                self.list.SetItem(citem)
                self.list.SetStringItem(i, 1, item[0])
            else:
                self.list.InsertStringItem(i, item[0])
            for j,item in enumerate(item[1:len(colLabels)]):
                item = wordwrap(StripIndents(item,True), colWidths[j+1], wx.ClientDC(self))
                item += "\n==========================================="
                self.list.SetStringItem(i,1+j+inc, item)
        # make buttons
        mainSizer.Add(self.list, 1, wx.EXPAND|wx.ALL, 1)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        OKbtn = wx.Button(self, wx.ID_OK)
        OKbtn.SetDefault()
        btnsizer.Add(OKbtn)
        if not checkLbl:
            btn = wx.Button(self, wx.ID_CLOSE,"Cancel") 
            btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        # bindings for close of window, double-click,...
        self.Bind(wx.EVT_CLOSE, self._onClose)
        if not checkLbl:
            OKbtn.Bind(wx.EVT_BUTTON,self._onSelect)
            self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self._onSelect)
            btn.Bind(wx.EVT_BUTTON,self._onClose)            
        self.SetSizer(mainSizer)
        self.ShowModal()
    def _onClose(self,event):
        event.Skip()
        self.EndModal(wx.ID_CANCEL)
    def _onSelect(self,event):
        if self.list.GetNextSelected(-1) == -1: return
        self.Selection = self.list.GetNextSelected(-1)
        self.EndModal(wx.ID_OK)

def MultiColMultiSelDlg(parent, title, header, colInfo, choices):
    '''Provides a dialog widget that can be used to select multiple items
    from a multicolumn list. 
    
    :param wx.Frame parent: the parent frame (or None)
    :param str title: A title for the dialog window
    :param str header: A instruction string for the dialog window
    :param list colInfo: contains three items for each column: a label for the column, 
      a width for the column (in pixels), and True if the column should be right justified.
    :param list choices: a nested list with values for each row in the table. Within each row
      should be a list of values for each column. There must be at least one value, but it is
      OK to have more or fewer values than there are column labels (colInfo). Extra are ignored
      and unspecified columns are left blank.
    :returns: a list of bool values for each entry in choices, True if selected, or
      None is the dialog is cancelled.
    
    Example use::

      choices = [('xmltodict', 'Bruker .brml Importer'),
                 ('zarr', 'MIDAS Zarr importer'),
                 ('h5py', 'HDF5 image importer'),
                 ('hdf5', 'HDF5 image importer')]
      colInfo = [('package', 50, False),
                 ('needed by', 200, True)]
      res = G2G.MultiColMultiSelDlg(parent, 'window title', 'Instructions', colInfo, choices)
    '''
    dlg = wx.Dialog(parent,wx.ID_ANY,title,
        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    txt = wx.StaticText(dlg,wx.ID_ANY,header)
    txt.Wrap(300)
    mainSizer.Add(txt)
    lst = wx.ListCtrl(dlg, wx.ID_ANY, style=wx.LC_REPORT)
    lst.EnableCheckBoxes()
    lst.InsertColumn(0, 'Sel')
    lst.SetColumnWidth(0, 30)
    cols = len(colInfo)
    for i,(lbl,wid,rgt) in enumerate(colInfo):
        if rgt:
            lst.InsertColumn(i+1, lbl, wx.LIST_FORMAT_RIGHT)
        else:
            lst.InsertColumn(i+1, lbl)
        if type(wid) is int:
            lst.SetColumnWidth(i+1, wid)
        else:
            lst.SetColumnWidth(i+1, wx.LIST_AUTOSIZE)
    for line in choices:
        index = lst.InsertItem(lst.GetItemCount(),'')
        for i,lbl in enumerate(line[:cols]):
            lst.SetItem(index, i+1, lbl)
    mainSizer.Add(lst,1,wx.EXPAND,1)
    btnsizer = wx.StdDialogButtonSizer()
    btn = wx.Button(dlg, wx.ID_OK, 'Install Selected')
    btn.SetDefault()
    btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_OK))
    btnsizer.AddButton(btn)
    btn = wx.Button(dlg, wx.ID_CANCEL)
    btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_CANCEL))
    btnsizer.AddButton(btn)
    btnsizer.Realize()
    mainSizer.Add(btnsizer, 0, wx.EXPAND|wx.ALL, 5)
    dlg.SetSizer(mainSizer)
    dlg.CenterOnParent()
    try:
        if dlg.ShowModal() == wx.ID_OK:
            return [lst.IsItemChecked(i) for i,c in enumerate(choices)]
        return
    finally:
        dlg.Destroy()
        
################################################################################
class OrderBox(wxscroll.ScrolledPanel):
    '''Creates a panel with scrollbars where items can be ordered into columns
    
    :param list keylist: is a list of keys for column assignments
    :param dict vallookup: is a dict keyed by names in keylist where each item is a dict. 
      Each inner dict contains variable names as keys and their associated values 
    :param dict posdict: is a dict keyed by names in keylist where each item is a dict. 
      Each inner dict contains column numbers as keys and their associated
      variable name as a value. This is used for both input and output.
      
    '''
    def __init__(self,parent,keylist,vallookup,posdict,*arg,**kw):
        self.keylist = keylist
        self.vallookup = vallookup
        self.posdict = posdict
        self.maxcol = 0
        for nam in keylist:
            posdict = self.posdict[nam]
            if posdict.keys():
                self.maxcol = max(self.maxcol, max(posdict))
        wxscroll.ScrolledPanel.__init__(self,parent,wx.ID_ANY,*arg,**kw)
        self.GBsizer = wx.GridBagSizer(4,4)
        self.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
        self.SetSizer(self.GBsizer)
        colList = [str(i) for i in range(self.maxcol+2)]
        for i in range(self.maxcol+1):
            wid = wx.StaticText(self,wx.ID_ANY,str(i),style=wx.ALIGN_CENTER)
            wid.SetBackgroundColour(DULL_YELLOW)
            wid.SetMinSize((50,-1))
            self.GBsizer.Add(wid,(0,i),flag=wx.EXPAND)
        self.chceDict = {}
        for row,nam in enumerate(self.keylist):
            posdict = self.posdict[nam]
            for col in posdict:
                lbl = posdict[col]
                pnl = wx.Panel(self,wx.ID_ANY)
                pnl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE))
                insize = wx.BoxSizer(wx.VERTICAL)
                wid = wx.Choice(pnl,wx.ID_ANY,choices=colList)
                insize.Add(wid,0,wx.EXPAND|wx.BOTTOM,3)
                wid.SetSelection(col)
                self.chceDict[wid] = (row,col)
                wid.Bind(wx.EVT_CHOICE,self.OnChoice)
                wid = wx.StaticText(pnl,wx.ID_ANY,lbl)
                insize.Add(wid,0,flag=wx.EXPAND)
                try:
                    val = G2fil.FormatSigFigs(self.vallookup[nam][lbl],maxdigits=8)
                except KeyError:
                    val = '?'
                wid = wx.StaticText(pnl,wx.ID_ANY,'('+val+')')
                insize.Add(wid,0,flag=wx.EXPAND)
                pnl.SetSizer(insize)
                self.GBsizer.Add(pnl,(row+1,col),flag=wx.EXPAND)
        self.SetAutoLayout(1)
        self.SetupScrolling()
        self.SetMinSize((
            min(700,self.GBsizer.GetSize()[0]),
            self.GBsizer.GetSize()[1]+20))
    def OnChoice(self,event):
        '''Called when a column is assigned to a variable
        '''
        row,col = self.chceDict[event.EventObject] # which variable was this?
        newcol = event.Selection # where will it be moved?
        if newcol == col:
            return # no change: nothing to do!
        prevmaxcol = self.maxcol # save current table size
        key = self.keylist[row] # get the key for the current row
        lbl = self.posdict[key][col] # selected variable name
        lbl1 = self.posdict[key].get(col+1,'') # next variable name, if any
        # if a posXXX variable is selected, and the next variable is posXXX, move them together
        repeat = 1
        if lbl[:3] == 'pos' and lbl1[:3] == 'int' and lbl[3:] == lbl1[3:]:
            repeat = 2
        for i in range(repeat): # process the posXXX and then the intXXX (or a single variable)
            col += i
            newcol += i
            if newcol in self.posdict[key]:
                # find first non-blank after newcol
                for mtcol in range(newcol+1,self.maxcol+2):
                    if mtcol not in self.posdict[key]: break
                l1 = range(mtcol,newcol,-1)+[newcol]
                l = range(mtcol-1,newcol-1,-1)+[col]
            else:
                l1 = [newcol]
                l = [col]
            # move all of the items, starting from the last column
            for newcol,col in zip(l1,l):
                #print 'moving',col,'to',newcol
                self.posdict[key][newcol] = self.posdict[key][col]
                del self.posdict[key][col]
                self.maxcol = max(self.maxcol,newcol)
                obj = self.GBsizer.FindItemAtPosition((row+1,col))
                self.GBsizer.SetItemPosition(obj.GetWindow(),(row+1,newcol))
                for wid in obj.GetWindow().Children:
                    if wid in self.chceDict:
                        self.chceDict[wid] = (row,newcol)
                        wid.SetSelection(self.chceDict[wid][1])
        # has the table gotten larger? If so we need new column heading(s)
        if prevmaxcol != self.maxcol:
            for i in range(prevmaxcol+1,self.maxcol+1):
                wid = wx.StaticText(self,wx.ID_ANY,str(i),style=wx.ALIGN_CENTER)
                wid.SetBackgroundColour(DULL_YELLOW)
                wid.SetMinSize((50,-1))
                self.GBsizer.Add(wid,(0,i),flag=wx.EXPAND)
            colList = [str(i) for i in range(self.maxcol+2)]
            for wid in self.chceDict:
                wid.SetItems(colList)
                wid.SetSelection(self.chceDict[wid][1])
        self.GBsizer.Layout()
        self.FitInside()
        
################################################################################
def GetImportFile(G2frame, message, defaultDir="", defaultFile="",
    style=wx.FD_OPEN, parent=None,*args, **kwargs):
    '''Uses a customized dialog that gets files from the appropriate import directory. 
    Arguments are used the same as in :func:`wx.FileDialog`. Selection of
    multiple files is allowed if argument style includes wx.FD_MULTIPLE.

    The default initial directory (unless overridden with argument defaultDir)
    is found in G2frame.TutorialImportDir, config setting Import_directory or
    G2frame.LastImportDir, see :func:`GetImportPath`.

    The path of the first file entered is used to set G2frame.LastImportDir
    and optionally config setting Import_directory.

    :returns: a list of files or an empty list
    '''
    if not parent: parent = G2frame
    pth = GetImportPath(G2frame)
    #if GSASIIpath.GetConfigValue('debug'):
    #    print('debug: GetImportFile from '+defaultDir)
    #    print('debug: GetImportFile pth '+pth)
    dlg = wx.FileDialog(parent, message, defaultDir, defaultFile, *args,style=style, **kwargs)
#    dlg.CenterOnParent()
    if not defaultDir and pth: dlg.SetDirectory(pth)
    try:
        if dlg.ShowModal() == wx.ID_OK:
            if style & wx.FD_MULTIPLE:
                filelist = dlg.GetPaths()
                if len(filelist) == 0: return []
            else:
                filelist = [dlg.GetPath(),]
            # not sure if we want to do this (why use wx.CHANGE_DIR?)
            if style & wx.FD_CHANGE_DIR: # to get Mac/Linux to change directory like windows!
                os.chdir(dlg.GetDirectory())
        else: # cancel was pressed
            return []
    finally:
        dlg.Destroy()
    # save the path of the first file and reset the TutorialImportDir variable
    pth = os.path.split(os.path.abspath(filelist[0]))[0]
    if GSASIIpath.GetConfigValue('Save_paths'): SaveImportDirectory(pth)
    G2frame.LastImportDir = pth
    G2frame.TutorialImportDir = None
    return filelist

def GetImportPath(G2frame):
    '''Determines the default location to use for importing files. Tries sequentially
    G2frame.TutorialImportDir, config var Import_directory, G2frame.LastImportDir
    and G2frame.LastGPXdir
    
    :returns: a string containing the path to be used when reading files or '.'
      if none of the above are specified.
    '''
    if G2frame.TutorialImportDir:
        if os.path.exists(G2frame.TutorialImportDir):
            return G2frame.TutorialImportDir
        elif GSASIIpath.GetConfigValue('debug'):
            print('DBG_Tutorial location (TutorialImportDir) not found: '+G2frame.TutorialImportDir)
    pth = GSASIIpath.GetConfigValue('Import_directory')
    if pth:
        pth = os.path.expanduser(pth)
        if os.path.exists(pth):
            return pth
        elif GSASIIpath.GetConfigValue('debug'):
            print('Ignoring Config Import_directory value: '+GSASIIpath.GetConfigValue('Import_directory'))
    if G2frame.LastImportDir:
        if os.path.exists(G2frame.LastImportDir):
            return G2frame.LastImportDir
        elif GSASIIpath.GetConfigValue('debug'):
            print('DBG_Warning: G2frame.LastImportDir not found = '+G2frame.LastImportDir)
    elif G2frame.LastGPXdir:
        return G2frame.LastGPXdir
    print('Import path not found - set to current directory')      #now shouldn't happen
    return '.'

def GetExportPath(G2frame):
    '''Determines the default location to use for writing files. Tries sequentially
    G2frame.LastExportDir and G2frame.LastGPXdir.
    
    :returns: a string containing the path to be used when writing files or '.'
      if none of the above are specified.
    '''
    if G2frame.LastExportDir:
        return G2frame.LastExportDir
    elif G2frame.LastGPXdir:
        return G2frame.LastGPXdir
    print('Export path not found - set to current directory')      #now shouldn't happen
    return '.'


################################################################################
class SGMessageBox(wx.Dialog):
    ''' Special version of MessageBox that displays space group & super space group text
    in two blocks
    '''
    def __init__(self,parent,title,text,table,spins=[],):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title,pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.text = text
        self.table = table
        self.panel = wx.Panel(self)
        self.spins = spins
        self.useAlt = False
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((0,10))
        for line in text:
            mainSizer.Add(wx.StaticText(self.panel,label='     %s     '%(line)))
        ncol = self.table[0].count(',')+1
        tableSizer = wx.FlexGridSizer(0,2*ncol+3,0,0)
        j = 0
        for item in self.table:
            if 'for' in item:
                mainSizer.Add(tableSizer,0,wx.ALIGN_LEFT)
                mainSizer.Add(wx.StaticText(self.panel,label=item),0)
                tableSizer = wx.FlexGridSizer(0,2*ncol+3,0,0)
                continue
            num,flds = item.split(')')
            tableSizer.Add(wx.StaticText(self.panel,label='     %s  '%(num+')')),0,WACV|wx.ALIGN_LEFT)            
            flds = flds.replace(' ','').split(',')
            for i,fld in enumerate(flds):
                if i < ncol-1:
                    text = wx.StaticText(self.panel,label='%s, '%(fld))
                else:
                    text = wx.StaticText(self.panel,label='%s'%(fld))
                if len(self.spins) and self.spins[j] < 0:
                    text.SetForegroundColour('Red')
                tableSizer.Add(text,0,WACV|wx.ALIGN_RIGHT)
            if not j%2:
                tableSizer.Add((20,0))
            j += 1
            
        def OnPrintOps(event):
            print(' Symmetry operations for %s:'%self.text[0].split(':')[1])
            for iop,opText in enumerate(G2spc.TextOps(self.text,self.table,reverse=True)):
                if self.useAlt:
                    opText = opText.replace('x','x1').replace('y','x2').replace('z','x3').replace('t','x4')
                if len(self.spins):
                    if self.spins[iop] > 0:
                        print('%s,+%d'%(opText.replace(' ',''),self.spins[iop]))
                    else:
                        print('%s,%d'%(opText.replace(' ',''),self.spins[iop]))
                else:    
                    print(opText.replace(' ',''))
                    
        def OnAlt(event):
            self.useAlt = altBtn.GetValue()
            
        mainSizer.Add(tableSizer,0,wx.ALIGN_LEFT)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        OKbtn = wx.Button(self.panel, wx.ID_OK)
        OKbtn.Bind(wx.EVT_BUTTON, self.OnOk)
        btnsizer.Add(OKbtn)
        printBtn = wx.Button(self.panel,label='Print Ops')
        printBtn.Bind(wx.EVT_BUTTON, OnPrintOps)
        btnsizer.Add(printBtn)
        altBtn = wx.CheckBox(self.panel,label=' Use alt. symbols?')
        altBtn.Bind(wx.EVT_CHECKBOX,OnAlt)
        btnsizer.Add(altBtn,0,WACV)
        mainSizer.Add((0,10))
        mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
        size = self.GetSize()
        self.SetSize([size[0]+20,size[1]])

    def Show(self):
        '''Use this method after creating the dialog to post it
        '''
        self.ShowModal()
        return

    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)

################################################################################
class SGMagSpinBox(wx.Dialog):
    ''' Special version of MessageBox that displays magnetic spin text
    '''
    def __init__(self,parent,title,text,table,Cents,names,spins,ifGray):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title,pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER,size=wx.Size(420,350))
        self.text = text
        self.table = table
        self.names = names
        Nnames = len(self.names)
        self.spins = spins
        self.ifGray = ifGray
        self.PrintTable = [' Magnetic symmetry operations for %s:'%self.text[0].split(':')[1],]
        self.panel = wxscroll.ScrolledPanel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((0,10))
        cents = [0,]
        if len(Cents) > 1:
            cents = self.text[-1].split(';')
        for line in self.text:
            mainSizer.Add(wx.StaticText(self.panel,label='     %s     '%(line)),0)
            if 'equivalent' in line:
                break
        ncol = self.table[0].count(',')+2
        nG = 1
        j = 0
        for ng in range(nG):
            if ng:
                mainSizer.Add(wx.StaticText(self.panel,label="      for (0,0,0)+1'"),0)
                j = 0
            for ic,cent in enumerate(cents):
                Cent = np.zeros(3)
                if cent:
                    cent = cent.strip(' (').strip(')+\n')
                    Cent = np.array(eval(cent)[:3])
#                Cent = np.array(Cents[ic])
                if ic:
                    if cent: cent = cent.strip(' (').strip(')+\n')
                    label = '      for (%s)+'%(cent)
                    if ng:     #test for gray operators
                        label += "1'"
                    mainSizer.Add(wx.StaticText(self.panel,label=label),0)
                tableSizer = wx.FlexGridSizer(0,2*ncol+3,0,0)
                for item in self.table:
                    if ')' not in item:
                        continue
                    flds = item.split(')')[1]
                    tableSizer.Add(wx.StaticText(self.panel,label='  (%2d)  '%(j+1)),0,WACV)            
                    flds = flds.replace(' ','').split(',')
                    for i,fld in enumerate(flds):
                        if i < ncol-1:
                            text = wx.StaticText(self.panel,label='%s, '%(fld))
                        else:
                            text = wx.StaticText(self.panel,label='%s '%(fld))
                        tableSizer.Add(text,0,WACV)
                    text = wx.StaticText(self.panel,label=' (%s) '%(self.names[j%Nnames]))
                    try:
                        if self.spins[j] < 0:
                            text.SetForegroundColour('Red')
                            item += ',-1'
                        else:
                            item += ',+1'
                    except IndexError:
                        print(self.spins,j,self.names[j%Nnames])
                        item += ',+1'
                    M,T,S = G2spc.MagText2MTS(item.split(')')[1].replace(' ',''),CIF=False)
                    T = (T+Cent)%1.
                    item = G2spc.MT2text([M,T],reverse=True)
                    if S > 0:
                        item += ',+1'
                    else:
                        item += ',-1'
                    self.PrintTable.append(item.replace(' ','').lower())
                    tableSizer.Add(text,0,WACV)
                    if not j%2:
                        tableSizer.Add((20,0))
                    j += 1
                mainSizer.Add(tableSizer,0)
            
            
        def OnPrintOps(event):
            for item in self.PrintTable:
                print(item)
            
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        OKbtn = wx.Button(self.panel, wx.ID_OK)
        btnsizer.Add(OKbtn)
        printBtn = wx.Button(self.panel,label='Print Ops')
        printBtn.Bind(wx.EVT_BUTTON, OnPrintOps)
        btnsizer.Add(printBtn)
        OKbtn.SetFocus()
        mainSizer.Add((0,10))
        mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER)
        
        self.panel.SetSizer(mainSizer)
        self.panel.SetAutoLayout(True)
        self.panel.SetScrollRate(10,10)
        self.panel.SendSizeEvent()


    def Show(self):
        '''Use this method after creating the dialog to post it
        '''
        self.ShowModal()
        return
    

################################################################################
class DisAglDialog(wx.Dialog):
    '''Distance/Angle Controls input dialog. After
    :meth:`ShowModal` returns, the results are found in
    dict :attr:`self.data`, which is accessed using :meth:`GetData`.

    :param wx.Frame parent: reference to parent frame (or None)
    :param dict data: a dict containing the current
      search ranges or an empty dict, which causes default values
      to be used.
      Will be used to set element `DisAglCtls` in 
      :ref:`Phase Tree Item <Phase_table>`
    :param dict default:  A dict containing the default
      search ranges for each element.
    :param bool Reset: if True (default), show Reset button
    :param bool Angle: if True (default), show angle radii
    '''
    def __init__(self,parent,data,default,Reset=True,Angle=True):
        text = 'Distance Angle Controls'
        if not Angle:
            text = 'Distance Controls'
        wx.Dialog.__init__(self,parent,wx.ID_ANY,text, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.default = default
        self.Reset = Reset
        self.Angle = Angle
        self.panel = None
        self._default(data,self.default)
        self.Draw(self.data)
        self.CenterOnParent()
                
    def _default(self,data,default):
        '''Set starting values for the search values, either from
        the input array or from defaults, if input is null
        '''
        if data:
            self.data = copy.deepcopy(data) # don't mess with originals
        else:
            self.data = {}
            self.data['Name'] = default['Name']
            self.data['Factors'] = [0.85,0.85]
            self.data['AtomTypes'] = default['AtomTypes']
            self.data['BondRadii'] = default['BondRadii'][:]
            self.data['AngleRadii'] = default['AngleRadii'][:]

    def Draw(self,data):
        '''Creates the contents of the dialog. Normally called
        by :meth:`__init__`.
        '''
        if self.panel: self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,'Controls for phase '+data['Name']),0,wx.LEFT,10)
        mainSizer.Add((10,10),1)
        
        ncol = 3
        if not self.Angle:
            ncol=2
        radiiSizer = wx.FlexGridSizer(0,ncol,5,5)
        radiiSizer.Add(wx.StaticText(self.panel,-1,' Type'),0,WACV)
        radiiSizer.Add(wx.StaticText(self.panel,-1,'Bond radii'),0,WACV)
        if self.Angle:
            radiiSizer.Add(wx.StaticText(self.panel,-1,'Angle radii'),0,WACV)
        self.objList = {}
        for id,item in enumerate(self.data['AtomTypes']):
            radiiSizer.Add(wx.StaticText(self.panel,-1,' '+item),0,WACV)
            bRadii = ValidatedTxtCtrl(self.panel,data['BondRadii'],id,nDig=(10,3))
            radiiSizer.Add(bRadii,0,WACV)
            if self.Angle:
                aRadii = ValidatedTxtCtrl(self.panel,data['AngleRadii'],id,nDig=(10,3))
                radiiSizer.Add(aRadii,0,WACV)
        mainSizer.Add(radiiSizer,0,wx.EXPAND)
        Names = ['Bond']
        factorSizer = wx.FlexGridSizer(0,2,5,5)
        if self.Angle:
            Names = ['Bond','Angle']
        for i,name in enumerate(Names):
            factorSizer.Add(wx.StaticText(self.panel,-1,name+' search factor'),0,WACV)
            bondFact = ValidatedTxtCtrl(self.panel,data['Factors'],i,nDig=(10,3))
            factorSizer.Add(bondFact)
        mainSizer.Add(factorSizer,0,wx.EXPAND)
        
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        if self.Reset:
            ResetBtn = wx.Button(self.panel,-1,'Reset')
            ResetBtn.Bind(wx.EVT_BUTTON, self.OnReset)
            btnSizer.Add(ResetBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
    
    def GetData(self):
        'Returns the values from the dialog'
        return self.data
        
    def OnOk(self,event):
        'Called when the OK button is pressed'
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnReset(self,event):
        'Called when the Reset button is pressed'
        data = {}
        self._default(data,self.default)
        wx.CallAfter(self.Draw,self.data)
                
################################################################################
class ShowLSParms(wx.Dialog):
    '''Create frame to show least-squares parameters
    '''
    def __init__(self,G2frame,title,parmDict,varyList,fullVaryList,
                     Controls, size=(650,430)):
        
        wx.Dialog.__init__(self,G2frame,wx.ID_ANY,title,size=size,
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.parmChoice = 'Phase'
        self.G2frame = G2frame
        self.parmDict = parmDict
        self.varyList = varyList
        self.fullVaryList = fullVaryList
        self.Controls = Controls
        self.choiceDict = {}

        parmFrozen = Controls.get('parmFrozen',{})
        if G2frame.testSeqRefineMode():
            frozenList = set()
            for h in parmFrozen:
                if h == 'FrozenList': continue
                frozenList = frozenList.union(parmFrozen[h])
            self.frozenList = list(frozenList)
        elif 'FrozenList' in parmFrozen:
            self.frozenList = copy.copy(parmFrozen['FrozenList'])
        else:
            self.frozenList = []
        # make lists of variables of different types along with lists of parameter names, histogram #s, phase #s,...
        self.parmNames = sorted(list(parmDict.keys()))
        if '2' in platform.python_version_tuple()[0]: 
            basestr = basestring
        else:
            basestr = str
        splitNames = [item.split(':') for item in self.parmNames if len(item) > 3 and not isinstance(self.parmDict[item],basestr)]
        globNames = [':'.join(item) for item in splitNames if not item[0] and not item[1]]
        if len(globNames):
            self.choiceDict['Global'] = G2obj.SortVariables(globNames)
        self.globVars = sorted(list(set([' ',]+[item[2] for item in splitNames if not item[0] and not item[1]])))
        hisNames = [':'.join(item) for item in splitNames if not item[0] and item[1]]
        self.choiceDict['Histogram'] = G2obj.SortVariables(hisNames)
        self.hisNums = sorted(list(set([int(item.split(':')[1]) for item in hisNames])))
        self.hisNums = ['*',]+[str(item) for item in self.hisNums]
        self.hisVars = sorted(list(set([' ',]+[item[2] for item in splitNames if not item[0]])))
        phasNames = [':'.join(item) for item in splitNames if not item[1] and not item[2].startswith('is')]
        self.choiceDict['Phase'] = G2obj.SortVariables(phasNames)
        self.phasNums = sorted(['*',]+list(set([item.split(':')[0] for item in phasNames])))
        if '' in self.phasNums: self.phasNums.remove('')
        self.phasVars = sorted(list(set([' ',]+[item[2] for item in splitNames if not item[1] and not item[2].startswith('is')])))
        hapNames = [':'.join(item) for item in splitNames if item[0] and item[1]]
        self.choiceDict['Phase/Histo'] = G2obj.SortVariables(hapNames)
        self.hapVars = sorted(list(set([' ',]+[item[2] for item in splitNames if item[0] and item[1]])))
        
        self.hisNum = '*'
        self.phasNum = '*'
        self.varName = ' '
        self.listSel = 'Refined'
        self.DrawPanel()
        
    def repaintScrollTbl(self):
        '''Shows the selected variables in a ListCtrl
        '''
        #start = time.time()
        self.varBox.SetContents(self)
        self.SendSizeEvent()
        #if GSASIIpath.GetConfigValue('debug'):
        #    print('repaintScrollTbl',time.time()-start)
                    
    def DrawPanel(self):
        '''Draws the contents of the entire dialog. Called initially & when radio buttons are pressed
        '''
        def _OnParmSel(event):
            'New parameter type, reset var name choice as list changes'
            self.parmChoice = parmSel.GetStringSelection()
            if varSel:
                varSel.SetSelection(0)
                self.varName = ' '
            wx.CallLater(100,self.DrawPanel)
            
        def OnPhasSel(event):
            'phase has been selected'
            event.Skip()
            self.phasNum = phasSel.GetValue()
            if varSel:
                try:
                    varSel.SetSelection(varSel.GetItems().index(self.varName))
                except:
                    varSel.SetSelection(0)
                    self.varName = ' '
            wx.CallAfter(self.repaintScrollTbl)

        def OnHistSel(event):
            'histogram has been selected'
            event.Skip()
            self.hisNum = histSel.GetValue()
            if varSel:
                try:
                    varSel.SetSelection(varSel.GetItems().index(self.varName))
                except:
                    varSel.SetSelection(0)
                    self.varName = ' '
            wx.CallAfter(self.repaintScrollTbl)
            
        def OnVarSel(event):
            'parameter name has been selected'
            event.Skip()
            self.varName = varSel.GetValue()
            if phasSel:
                try:
                    phasSel.SetSelection(phasSel.GetItems().index(self.phasNum))
                except:
                    phasSel.SetSelection(0)
                    self.phasNum = '*'
            if histSel:
                try:
                    histSel.SetSelection(histSel.GetItems().index(self.hisNum))
                except:
                    histSel.SetSelection(0)
                    self.hisNum = '*'
            wx.CallAfter(self.repaintScrollTbl)
            
        def OnListSel(event):
            self.listSel = listSel.GetStringSelection()
            wx.CallLater(100,self.DrawPanel)
                        
        def OnVarSpin(event):
            '''Respond when any of the SpinButton widgets are pressed'''
            event.Skip()
            Spinner = event.GetEventObject()
            move = Spinner.GetValue()
            Spinner.SetValue(0)
            varSel,binding = self.SpinDict[Spinner.GetId()]
            i = varSel.GetSelection() - move
            if i < 0:
                i = varSel.GetCount()-1
            elif i >= varSel.GetCount():
                i = 0
            varSel.SetSelection(i)
            wx.CallLater(100,binding,event)

        def AddSpinner(varSizer,label,SelCtrl,binding):
            '''Add a label and a SpinButton to a Combo widget (SelCtrl)
            Saves a pointer to the combo widget and the callback used by that widget
            '''
            SelCtrl.Bind(wx.EVT_COMBOBOX,binding)
            varSizer.Add(wx.StaticText(self,label=label))
            varSelSizer = wx.BoxSizer(wx.HORIZONTAL)
            varSelSizer.Add(SelCtrl,0)
            varSpin = wx.SpinButton(self,style=wx.SP_VERTICAL)
            varSpin.SetValue(0)
            varSpin.SetRange(-1,1)
            varSpin.Bind(wx.EVT_SPIN, OnVarSpin)
            self.SpinDict[varSpin.GetId()] = SelCtrl,binding
            varSelSizer.Add(varSpin,0)
            varSizer.Add(varSelSizer,0)

        if self.GetSizer(): self.GetSizer().Clear(True)
        self.SpinDict = {}
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        num = len(self.varyList)
        mainSizer.Add(wx.StaticText(self,label='View Parameters in Project'),0,wx.ALIGN_CENTER)
        parmSizer = wx.BoxSizer(wx.HORIZONTAL)
        parmSizer.Add(wx.StaticText(self,label=' Number of refined variables: {}'.format(num)),0,wx.ALIGN_LEFT)
        if len(self.varyList) != len(self.fullVaryList):
            num = len(self.fullVaryList) - len(self.varyList)
            parmSizer.Add(wx.StaticText(self,label=
                ' + {} varied via constraints'.format(
                    len(self.fullVaryList) - len(self.varyList))
                                        ))
        parmFrozen = self.Controls.get('parmFrozen',{})
        fcount = 0
        if self.G2frame.testSeqRefineMode():
            for h in parmFrozen:
                if h == 'FrozenList': continue
                fcount += len(parmFrozen[h])
        elif 'FrozenList' in parmFrozen:
            fcount = len(parmFrozen['FrozenList'])
        if fcount:
            parmSizer.Add(wx.StaticText(self,label=
                ' - {} frozen variables'.format(fcount)))
        mainSizer.Add(parmSizer)
        choice = ['Phase','Phase/Histo','Histogram']
        if 'Global' in self.choiceDict:
            choice += ['Global',]
        parmSizer = wx.BoxSizer(wx.HORIZONTAL)
        parmSel = wx.RadioBox(self,wx.ID_ANY,'Parameter type:',choices=choice,
            majorDimension=1,style=wx.RA_SPECIFY_COLS)
        parmSel.Bind(wx.EVT_RADIOBOX,_OnParmSel)
        parmSel.SetStringSelection(self.parmChoice)
        parmSizer.Add(parmSel,0)
        
        selectionsSizer = wx.BoxSizer(wx.VERTICAL)        
        varSizer = wx.BoxSizer(wx.VERTICAL)
        varSel = None
        if self.parmChoice != 'Global': 
            if self.parmChoice in ['Phase',]:
                varSel = wx.ComboBox(self,choices=self.phasVars,value=self.varName,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
            elif self.parmChoice in ['Histogram',]:
                varSel = wx.ComboBox(self,choices=self.hisVars,value=self.varName,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
            elif self.parmChoice in ['Phase/Histo',]:
                varSel = wx.ComboBox(self,choices=self.hapVars,value=self.varName,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
            AddSpinner(varSizer,'Parameter',varSel,OnVarSel)
        selectionsSizer.Add(varSizer,0)
                
        varSizer = wx.BoxSizer(wx.HORIZONTAL)
        phasSel = None
        if self.parmChoice in ['Phase','Phase/Histo'] and len(self.phasNums) > 1:
            numSizer = wx.BoxSizer(wx.VERTICAL)
            phasSel = wx.ComboBox(self,choices=self.phasNums,value=self.phasNum,
                style=wx.CB_READONLY|wx.CB_DROPDOWN,size=(50,-1))
            AddSpinner(numSizer,'Phase',phasSel,OnPhasSel)
            varSizer.Add(numSizer)

        histSel = None
        if self.parmChoice in ['Histogram','Phase/Histo'] and len(self.hisNums) > 1:
            numSizer = wx.BoxSizer(wx.VERTICAL)
            histSel = wx.ComboBox(self,choices=self.hisNums,value=self.hisNum,
                style=wx.CB_READONLY|wx.CB_DROPDOWN,size=(50,-1))
            AddSpinner(numSizer,'Histogram',histSel,OnHistSel)
            varSizer.Add(numSizer)
        selectionsSizer.Add(varSizer,0)
        parmSizer.Add(selectionsSizer,0)
        refChoices = ['All','Refined']
        txt = ('"R" indicates a refined variable\n'+
               '"C" indicates generated from a user entered constraint')
        if fcount:
            refChoices += ['Frozen']
            txt += '\n"F" indicates a variable that is Frozen due to exceeding min/max'
        
        listSel = wx.RadioBox(self,wx.ID_ANY,'Refinement Status:',
            choices=refChoices,
            majorDimension=0,style=wx.RA_SPECIFY_COLS)
        listSel.SetStringSelection(self.listSel)
        listSel.Bind(wx.EVT_RADIOBOX,OnListSel)
        parmSizer.Add(listSel,0,wx.CENTER|wx.ALL,15)
        mainSizer.Add(parmSizer,0)
        
        self.countSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(self.countSizer)
        self.headSizer = wx.BoxSizer(wx.HORIZONTAL) # non-scrolling header        
        mainSizer.Add(self.headSizer,0)
        self.varBox = VirtualVarBox(self)
        mainSizer.Add(self.varBox,1,wx.ALL|wx.EXPAND,1)
        mainSizer.Add(
            wx.StaticText(self,label=txt),0, wx.ALL,0)
        
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)          # make Close button 
        btn = wx.Button(self, wx.ID_CLOSE,"Close") 
        btn.Bind(wx.EVT_BUTTON,self._onClose)
        btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        self.SetSizer(mainSizer)
        wx.CallAfter(self.repaintScrollTbl)
                
    def _onClose(self,event):
        self.EndModal(wx.ID_CANCEL)

class VirtualVarBox(wx.ListCtrl):
    def __init__(self, parent):
        self.parmWin = parent
        #patch (added Oct 2020) convert variable names for parm limits to G2VarObj
        G2sc.patchControls(self.parmWin.Controls)
        # end patch
        wx.ListCtrl.__init__(
            self, parent, -1,
            style=wx.LC_REPORT|wx.LC_VIRTUAL|wx.LC_HRULES|wx.LC_VRULES
            )

        for i,(lbl,wid) in enumerate(zip(
                ('#', "Parameter", "Ref", "Value", "Min", "Max", "Explanation"),
                (40 , 125        , 30    ,  100   ,  75  ,  75  , 700),)):  
            self.InsertColumn(i, lbl)
            self.SetColumnWidth(i, wid)

        self.SetItemCount(0)

        try:
            self.attr1 = wx.ItemAttr()
        except:
            self.attr1 = wx.ListItemAttr() # deprecated in wx4.1
        self.attr1.SetBackgroundColour((255,255,150))

        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnRowSelected)

    def SetContents(self,parent):
        self.varList = []
        for name in parent.choiceDict[parent.parmChoice]:
            if '2' in platform.python_version_tuple()[0]: 
                basestr = basestring
            else:
                basestr = str
            if isinstance(parent.parmDict[name],basestr): continue
            if 'Refined' in parent.listSel and (name not in parent.fullVaryList
                                              ) and (name not in parent.varyList):
                continue
            if 'Frozen' in parent.listSel and not (
                    name in self.parmWin.fullVaryList and
                    name in self.parmWin.frozenList):
                continue
            if 'Phase' in parent.parmChoice:
                if parent.phasNum != '*' and name.split(':')[0] != parent.phasNum: continue
            if 'Histo' in parent.parmChoice:
                if parent.hisNum != '*' and name.split(':')[1] != parent.hisNum: continue
            if (parent.varName != ' ') and (parent.varName not in name): continue
            self.varList.append(name)
        oldlen = self.GetItemCount()
        self.SetItemCount(len(self.varList))
        
    def OnRowSelected(self, event, row=None):
        'Creates an edit window when a parameter is selected'
        def ResetFrozen(event):
            '''release a frozen parameter (from all histograms in the case of a 
            sequential fit).
            '''
            if name in self.parmWin.frozenList:
                del self.parmWin.frozenList[self.parmWin.frozenList.index(name)]
            parmFrozen = self.parmWin.Controls.get('parmFrozen',{})
            if self.parmWin.G2frame.testSeqRefineMode():
                for h in parmFrozen:
                    if h == 'FrozenList': continue
                    if name in parmFrozen[h]:
                        del parmFrozen[h][parmFrozen[h].index(name)]
            elif 'FrozenList' in parmFrozen:
                if name in parmFrozen['FrozenList']:
                    del parmFrozen['FrozenList'][parmFrozen['FrozenList'].index(name)]
            dlg.EndModal(wx.ID_CANCEL)
            self.parmWin.SendSizeEvent()

        def delM(event):
            'Get event info & prepare to delete Max or Min limit'
            if hasattr(event.EventObject,'max'):
                d = self.parmWin.Controls['parmMaxDict']
            else:
                d = self.parmWin.Controls['parmMinDict']
            # close sub-dialog then delete item and redraw
            dlg.EndModal(wx.ID_OK)
            wx.CallAfter(delMafter,d,name)
        def delMafter(d,name):
            'Delete Max or Min limit once dialog is deleted'
            key,val = G2obj.prmLookup(name,d) # is this a wild-card?
            if val is not None:
                del d[key]
            self.OnRowSelected(None, row)
            
        def AddM(event):
            'Get event info & add a Max or Min limit'
            if hasattr(event.EventObject,'max'):
                d = self.parmWin.Controls['parmMaxDict']
            else:
                d = self.parmWin.Controls['parmMinDict']
            # close sub-dialog then delete item and redraw
            dlg.EndModal(wx.ID_OK)
            wx.CallAfter(AddMafter,d,name)
        def AddMafter(d,name):
            'Add a Max or Min limit & redraw'
            try:
                d[G2obj.G2VarObj(name)] = float(value)
            except:
                pass
            self.OnRowSelected(None, row)
        def SetWild(event):
            'Get event info & prepare to set/clear item as wildcard'
            if hasattr(event.EventObject,'max'):
                d = self.parmWin.Controls['parmMaxDict']
            else:
                d = self.parmWin.Controls['parmMinDict']
            ns = name.split(':')
            if hasattr(event.EventObject,'hist'):
                ns[1] = '*'
            else:
                ns[3] = '*'
            wname = ':'.join(ns)            
            # close sub-dialog then delete item and redraw
            dlg.EndModal(wx.ID_OK)
            wx.CallAfter(SetWildAfter,d,name,wname,event.EventObject.GetValue())
        def SetWildAfter(d,name,wname,mode):
            'Set/clear item as wildcard & delete old name(s), redraw'
            n,val = G2obj.prmLookup(name,d) # is this a wild-card?
            if val is None:
                print('Error: Limit for parameter {} not found. Should not happen'.format(name))
                return
            if mode: # make wildcard
                for n in list(d.keys()): # delete names matching wildcard
                    if str(n) == wname: continue
                    if n == wname: # respects wildcards
                        del d[n]
                d[G2obj.G2VarObj(wname)] = val
            else:
                del d[n]
                d[G2obj.G2VarObj(name)] = val
            self.OnRowSelected(None, row)

        # start of OnRowSelected 
        if event is not None:
            row = event.Index
        elif row is None:
            print('Error: row and event should not both be None!')
            return
        name = self.varList[row]
        dlg = wx.Dialog(self.parmWin,wx.ID_ANY,'Parameter {} info'.format(name),
                            size=(600,-1),
                        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5))
        subSizer = wx.BoxSizer(wx.HORIZONTAL)
        subSizer.Add((-1,-1),1,wx.EXPAND)
        try:
            value = G2fil.FormatSigFigs(self.parmWin.parmDict[name])
        except TypeError:
            value = str(self.parmWin.parmDict[name])+' -?' # unexpected
        subSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
                        'Parameter "{}" information and settings. Value={}'
                                       .format(name,value)))
        subSizer.Add((-1,-1),1,wx.EXPAND)
        mainSizer.Add(subSizer,0,wx.EXPAND,0)
        mainSizer.Add((0,10))
        v = G2obj.getVarDescr(name)
        if v is not None and v[-1] is not None:
            txt = G2obj.fmtVarDescr(name)
            if txt:
                #txt = txt.replace('Ph=','Phase: ')
                #txt = txt.replace('Pwd=','Histogram: ')
                txtwid = wx.StaticText(dlg,wx.ID_ANY,'Parameter meaning is "'+txt+'"')
                txtwid.Wrap(580)
                mainSizer.Add(txtwid)
                mainSizer.Add((0,10))

        freezebtn = None
        if name in self.parmWin.fullVaryList and name in self.parmWin.frozenList:
            msg = "Parameter {} exceeded limits and has been frozen".format(name)
            freezebtn = wx.Button(dlg, wx.ID_ANY,'Unfreeze')
            freezebtn.Bind(wx.EVT_BUTTON, ResetFrozen)
        elif name in self.parmWin.varyList:
            msg = "Parameter {} is refined".format(name)
        elif name in self.parmWin.fullVaryList:
            msg = "Parameter {} is refined via a constraint".format(name)
        else:
            msg = ""
        if msg:
            subSizer = wx.BoxSizer(wx.HORIZONTAL)
            subSizer.Add(wx.StaticText(dlg,wx.ID_ANY,msg),0,wx.CENTER)
            if freezebtn:
                subSizer.Add(freezebtn,0,wx.ALL|wx.CENTER,5)
            mainSizer.Add(subSizer,0)

        # draw min value widgets
        mainSizer.Add((-1,10),0)
        if name not in self.parmWin.varyList and name in self.parmWin.fullVaryList:
            mainSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Limits not allowed on constrained variables'),0)
            for key in 'parmMinDict','parmMaxDict':
                d = self.parmWin.Controls[key]
                n,v = G2obj.prmLookup(name,d)
                if v is not None and str(n) == name:
                    try:  # strange hard to reproduce problem with this not working 
                        del d[n]
                    except:
                        if GSASIIpath.GetConfigValue('debug'):
                            print('debug: failed to delete ',name,'in',key)

        else:
            n,val = G2obj.prmLookup(name,self.parmWin.Controls['parmMinDict']) # is this a wild-card?
            if val is None: 
                addMbtn = wx.Button(dlg, wx.ID_ANY,'Add Lower limit')
                addMbtn.Bind(wx.EVT_BUTTON, AddM)
                mainSizer.Add(addMbtn,0)
            else:
                subSizer = wx.BoxSizer(wx.HORIZONTAL)
                subSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Minimum limit'),0,wx.CENTER)
                subSizer.Add(ValidatedTxtCtrl(dlg,self.parmWin.Controls['parmMinDict'],n,nDig=(10,2,'g')),0,WACV)
                delMbtn = wx.Button(dlg, wx.ID_ANY,'Delete',style=wx.BU_EXACTFIT)
                subSizer.Add((5,-1),0,WACV)
                subSizer.Add(delMbtn,0,WACV)
                delMbtn.Bind(wx.EVT_BUTTON, delM)
                if name.split(':')[1]:             # is this using a histogram?
                    subSizer.Add((5,-1),0,WACV)
                    wild = wx.CheckBox(dlg,wx.ID_ANY,label='Match all histograms ')
                    wild.SetValue(str(n).split(':')[1] == '*')
                    wild.Bind(wx.EVT_CHECKBOX,SetWild)
                    wild.hist = True
                    subSizer.Add(wild,0,WACV)
                elif len(name.split(':')) > 3:
                    subSizer.Add((5,-1),0,WACV)
                    wild = wx.CheckBox(dlg,wx.ID_ANY,label='Match all atoms ')
                    wild.SetValue(str(n).split(':')[3] == '*')
                    wild.Bind(wx.EVT_CHECKBOX,SetWild)
                    subSizer.Add(wild,0,WACV)
                mainSizer.Add(subSizer,0)
            # draw max value widgets
            mainSizer.Add((-1,10),0)
            n,val = G2obj.prmLookup(name,self.parmWin.Controls['parmMaxDict']) # is this a wild-card?
            if val is None: 
                addMbtn = wx.Button(dlg, wx.ID_ANY,'Add Upper limit')
                addMbtn.Bind(wx.EVT_BUTTON, AddM)
                addMbtn.max = True
                mainSizer.Add(addMbtn,0)
            else:
                subSizer = wx.BoxSizer(wx.HORIZONTAL)
                subSizer.Add(wx.StaticText(dlg,wx.ID_ANY,'Maximum limit'),0,wx.CENTER)
                subSizer.Add(ValidatedTxtCtrl(dlg,self.parmWin.Controls['parmMaxDict'],n,nDig=(10,2,'g')),0,WACV)
                delMbtn = wx.Button(dlg, wx.ID_ANY,'Delete',style=wx.BU_EXACTFIT)
                subSizer.Add((5,-1),0,WACV)
                subSizer.Add(delMbtn,0,WACV)
                delMbtn.Bind(wx.EVT_BUTTON, delM)
                delMbtn.max = True
                if name.split(':')[1]:             # is this using a histogram?
                    subSizer.Add((5,-1),0,WACV)
                    wild = wx.CheckBox(dlg,wx.ID_ANY,label='Match all histograms ')
                    wild.SetValue(str(n).split(':')[1] == '*')
                    wild.Bind(wx.EVT_CHECKBOX,SetWild)
                    wild.max = True
                    wild.hist = True
                    subSizer.Add(wild,0,WACV)
                elif len(name.split(':')) > 3:
                    subSizer.Add((5,-1),0,WACV)
                    wild = wx.CheckBox(dlg,wx.ID_ANY,label='Match all atoms ')
                    wild.SetValue(str(n).split(':')[3] == '*')
                    wild.Bind(wx.EVT_CHECKBOX,SetWild)
                    wild.max = True
                    subSizer.Add(wild,0,WACV)
                mainSizer.Add(subSizer,0)
            
        btnsizer = wx.StdDialogButtonSizer()
        OKbtn = wx.Button(dlg, wx.ID_OK)
        OKbtn.SetDefault()
        OKbtn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_OK))
        btnsizer.AddButton(OKbtn)
        btnsizer.Realize()
        mainSizer.Add((-1,5),1,wx.EXPAND,1)
        mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER,0)
        mainSizer.Add((-1,10))

        dlg.SetSizer(mainSizer)
        dlg.CenterOnParent()
        if dlg.ShowModal() != wx.ID_OK: # if not OK, destroy & reopen
            dlg.Destroy()
            wx.CallAfter(self.OnRowSelected, None, row)
            return
        dlg.Destroy()
        self.parmWin.SendSizeEvent()
        
    #-----------------------------------------------------------------
    # Callbacks to display info in table
    def OnGetItemText(self, item, col):
        name = self.varList[item]
        if col == 0:
            return str(item)
        elif col == 1:
            return name
        elif col == 2:
            if name in self.parmWin.fullVaryList and name in self.parmWin.frozenList:
                    return "F"
            elif name in self.parmWin.varyList:
                return "R"
            elif name in self.parmWin.fullVaryList:
                return "C"
            return ""
        elif col == 3:
            try:
                value = G2fil.FormatSigFigs(self.parmWin.parmDict[name])
            except TypeError:
                value = str(self.parmWin.parmDict[name])+' -?' # unexpected
            return value
        elif col == 4 or col == 5: # min/max value
            if col == 4: # min
                d = self.parmWin.Controls['parmMinDict']
            else:
                d = self.parmWin.Controls['parmMaxDict']
            n,val = G2obj.prmLookup(name,d)
            if val is None: return ""
            try:
                return G2fil.FormatSigFigs(val,8)
            except TypeError:
                return "?"
        elif col == 6:
            v = G2obj.getVarDescr(name)
            if v is not None and v[-1] is not None:
                txt = G2obj.fmtVarDescr(name)
                if txt: return txt
            return ""
        else:
            return "?"

    def OnGetItemAttr(self, item):
        name = self.varList[item]
        if name in self.parmWin.varyList and name in self.parmWin.frozenList:
            return self.attr1
        else:
            return None

#####  Customized Grid Support ################################################################################           
class GSGrid(wg.Grid):
    '''Basic wx.Grid implementation
    '''
    def __init__(self, parent, name=''):
        wg.Grid.__init__(self,parent,-1,name=name)
        if hasattr(parent.TopLevelParent,'currentGrids'):
            parent.TopLevelParent.currentGrids.append(self)      # save a reference to the grid in the Frame
        self.SetScrollRate(0,0)         #GSAS-II grids have no scroll bars by default
            
    def Clear(self):
        wg.Grid.ClearGrid(self)
        
    def SetCellReadOnly(self,r,c,readonly=True):
        self.SetReadOnly(r,c,isReadOnly=readonly)
        
    def SetCellStyle(self,r,c,color="white",readonly=True):
        self.SetCellBackgroundColour(r,c,color)
        self.SetReadOnly(r,c,isReadOnly=readonly)

    def SetTable(self, table, *args, **kwargs):
        '''Overrides the standard SetTable method with one that uses
        GridFractionEditor for all numeric columns (unless useFracEdit
        is false)
        '''
        setFracEdit = kwargs.get('useFracEdit',True)
        if 'useFracEdit' in kwargs: del kwargs['useFracEdit']
        wg.Grid.SetTable(self, table, *args, **kwargs)
        if setFracEdit:
            for i,t in enumerate(table.dataTypes):
                if not t.startswith(wg.GRID_VALUE_FLOAT): continue
                attr = wx.grid.GridCellAttr()
                attr.IncRef()
                attr.SetEditor(GridFractionEditor(self))
                self.SetColAttr(i, attr)

    def GetSelection(self):
        #this is to satisfy structure drawing stuff in G2plt when focus changes
        return None

    def InstallGridToolTip(self, rowcolhintcallback,
                           colLblCallback=None,rowLblCallback=None):
        '''code to display a tooltip for each item on a grid
        from http://wiki.wxpython.org/wxGrid%20ToolTips (buggy!), expanded to
        column and row labels using hints from
        https://groups.google.com/forum/#!topic/wxPython-users/bm8OARRVDCs

        :param function rowcolhintcallback: a routine that returns a text
          string depending on the selected row and column, to be used in
          explaining grid entries.
        :param function colLblCallback: a routine that returns a text
          string depending on the selected column, to be used in
          explaining grid columns (if None, the default), column labels
          do not get a tooltip.
        :param function rowLblCallback: a routine that returns a text
          string depending on the selected row, to be used in
          explaining grid rows (if None, the default), row labels
          do not get a tooltip.
        '''
        prev_rowcol = [None,None,None]
        def OnMouseMotion(event):
            # event.GetRow() and event.GetCol() would be nice to have here,
            # but as this is a mouse event, not a grid event, they are not
            # available and we need to compute them by hand.
            x, y = self.CalcUnscrolledPosition(event.GetPosition())
            row = self.YToRow(y)
            col = self.XToCol(x)
            hinttext = ''
            win = event.GetEventObject()
            if [row,col,win] == prev_rowcol: # no change from last position
                if event: event.Skip()
                return
            if win == self.GetGridWindow() and row >= 0 and col >= 0:
                hinttext = rowcolhintcallback(row, col)
            elif win == self.GetGridColLabelWindow() and col >= 0:
                if colLblCallback: hinttext = colLblCallback(col)
            elif win == self.GetGridRowLabelWindow() and row >= 0:
                if rowLblCallback: hinttext = rowLblCallback(row)
            else: # this should be the upper left corner, which is empty
                if event: event.Skip()
                return
            if hinttext is None: hinttext = ''
            if 'phoenix' in wx.version():
                win.SetToolTip(hinttext)
            else:
                win.SetToolTipString(hinttext)
            prev_rowcol[:] = [row,col,win]
            if event: event.Skip()
        if 'phoenix' in wx.version():
            self.GetGridWindow().Bind(wx.EVT_MOTION,OnMouseMotion)
            if colLblCallback: self.GetGridColLabelWindow().Bind(wx.EVT_MOTION,OnMouseMotion)
            if rowLblCallback: self.GetGridRowLabelWindow().Bind(wx.EVT_MOTION,OnMouseMotion)
        else:
            wx.EVT_MOTION(self.GetGridWindow(), OnMouseMotion)
            if colLblCallback: wx.EVT_MOTION(self.GetGridColLabelWindow(), OnMouseMotion)
            if rowLblCallback: wx.EVT_MOTION(self.GetGridRowLabelWindow(), OnMouseMotion)

    def setupPopup(self,lblList,callList):
        '''define a callback that creates a popup menu. The rows associated
        with the items selected items are selected in the table and if 
        an item is called from the menu, the corresponding function 
        is called to perform an action on the 

        :param list lblList: list of str items that will be placed in the 
          popup menu
        :param list callList: list of functions to be called when a 
        :returns: a callback that can be used to create the menu

        Sample usage::

            lblList = ('Delete','Set atom style','Set atom label',
                           'Set atom color','Set view point','Generate copy',
                           'Generate surrounding sphere','Transform atoms',
                           'Generate bonded')
            callList = (DrawAtomsDelete,DrawAtomStyle, DrawAtomLabel,
                            DrawAtomColor,SetViewPoint,AddSymEquiv,
                            AddSphere,TransformSymEquiv,
                            FillCoordSphere)
            onRightClick = drawAtoms.setupPopup(lblList,callList)
            drawAtoms.Bind(wg.EVT_GRID_CELL_RIGHT_CLICK, onRightClick)
            drawAtoms.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, onRightClick)

        '''
        def createPopup(event):
            def OnPopup(event):
                callback = callList[menuIndx.index(event.GetId())]
                self.ClearSelection()
                for r in indx:
                    self.SelectRow(r,True)
                callback(event)
            # get selections
            indx = self.GetSelectedRows()
            indx += [row for row,col in self.GetSelectedCells()]
            for top,bottom in zip([r for r,c in self.GetSelectionBlockTopLeft()],
                          [r for r,c in self.GetSelectionBlockBottomRight()]):
                indx += list(range(top,bottom+1))
            indx = list(set(indx))
            if len(indx) == 0: # nothing selected, get current row
                r,_ =  event.GetRow(),event.GetCol()
                if r < 0:
                    return
                indx = [r]
            # make a pop-up menu
            menu = wx.Menu()
            menuIndx = []
            for l in lblList:
                menuIndx.append(wx.NewIdRef())
                menu.Append(menuIndx[-1], l)
                self.Bind(wx.EVT_MENU, OnPopup, id=menuIndx[-1])
            self.PopupMenu(menu)
            menu.Destroy()
        return createPopup
        
    def completeEdits(self):
        'complete any outstanding edits'
        if self.IsCellEditControlEnabled(): # complete any grid edits in progress
            self.SaveEditControlValue()
            self.HideCellEditControl()
            self.DisableCellEditControl()
                
################################################################################           
class Table(wg.GridTableBase):
    '''Basic data table for use with GSgrid
    '''
    def __init__(self, data=[], rowLabels=None, colLabels=None, types = None):
        wg.GridTableBase.__init__(self)
        self.colLabels = colLabels
        self.rowLabels = rowLabels
        self.dataTypes = types
        self.data = data
        
    def AppendRows(self, numRows=1):
        self.data.append([])
        return True
        
    def CanGetValueAs(self, row, col, typeName):
        if self.dataTypes:
            colType = self.dataTypes[col].split(':')[0]
            if typeName == colType:
                return True
            else:
                return False
        else:
            return False

    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)

    def DeleteRow(self,pos):
        data = self.GetData()
        self.SetData([])
        new = []
        for irow,row in enumerate(data):
            if irow != pos:
                new.append(row)
        self.SetData(new)
        
    def GetColLabelValue(self, col):
        if self.colLabels:
            return self.colLabels[col]
            
    def GetData(self):
        data = []
        for row in range(self.GetNumberRows()):
            data.append(self.GetRowValues(row))
        return data
        
    def GetNumberCols(self):
        try:
            return len(self.colLabels)
        except TypeError:
            return None
        
    def GetNumberRows(self):
        return len(self.data)
        
    def GetRowLabelValue(self, row):
        if self.rowLabels:
            return self.rowLabels[row]
        
    def GetColValues(self, col):
        data = []
        for row in range(self.GetNumberRows()):
            data.append(self.GetValue(row, col))
        return data
        
    def GetRowValues(self, row):
        data = []
        for col in range(self.GetNumberCols()):
            data.append(self.GetValue(row, col))
        return data
        
    def GetTypeName(self, row, col):
        try:
            if self.data[row][col] is None:
                return wg.GRID_VALUE_STRING
            return self.dataTypes[col]
        except (TypeError,IndexError):
            return wg.GRID_VALUE_STRING

    def GetValue(self, row, col):
        try:
            if self.data[row][col] is None: return ""
            return self.data[row][col]
        except IndexError:
            return None
            
    def InsertRows(self, pos, rows):
        for row in range(rows):
            self.data.insert(pos,[])
            pos += 1
        
    def IsEmptyCell(self,row,col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
        
    def OnKeyPress(self, event):
        dellist = self.GetSelectedRows()
        if event.GetKeyCode() == wx.WXK_DELETE and dellist:
            grid = self.GetView()
            for i in dellist: grid.DeleteRow(i)
                
    def SetColLabelValue(self, col, label):
        numcols = self.GetNumberCols()
        if col > numcols-1:
            self.colLabels.append(label)
        else:
            self.colLabels[col]=label
        
    def SetData(self,data):
        for row in range(len(data)):
            self.SetRowValues(row,data[row])
                
    def SetRowLabelValue(self, row, label):
        self.rowLabels[row]=label
            
    def SetRowValues(self,row,data):
        self.data[row] = data
            
    def SetValue(self, row, col, value):
        def innerSetValue(row, col, value):
            try:
                self.data[row][col] = value
            except TypeError:
                return
            except IndexError: # has this been tested? 
                #print row,col,value
                if self.GetNumberRows() == 0: return
                # add a new row
                if row > self.GetNumberRows():
                    self.data.append([''] * self.GetNumberCols())
                elif col > self.GetNumberCols():
                    for row in range(self.GetNumberRows()): # bug fixed here
                        self.data[row].append('')
                #print self.data
                self.data[row][col] = value
        innerSetValue(row, col, value)

################################################################################
class GridFractionEditor(wg.PyGridCellEditor):
    '''A grid cell editor class that allows entry of values as fractions as well
    as sine and cosine values [as s() and c(), sin() or sind(), etc]. Any valid 
    Python expression will be evaluated. 

    The current value can be incremented, multiplied or divided by prefixing 
    an expression by +, * or / respectively. 
    '''
    def __init__(self,grid):
        if 'phoenix' in wx.version():
            wg.GridCellEditor.__init__(self)
        else:
            wg.PyGridCellEditor.__init__(self)

    def Create(self, parent, id, evtHandler):
        self._tc = wx.TextCtrl(parent, id, "")
        self._tc.SetInsertionPoint(0)
        self.SetControl(self._tc)

        if evtHandler:
            self._tc.PushEventHandler(evtHandler)

        self._tc.Bind(wx.EVT_CHAR, self.OnChar)

    def SetSize(self, rect):
        if 'phoenix' in wx.version():
            self._tc.SetSize(rect.x, rect.y, rect.width+2, rect.height+2,
                               wx.SIZE_ALLOW_MINUS_ONE)
        else:
            self._tc.SetDimensions(rect.x, rect.y, rect.width+2, rect.height+2,                                wx.SIZE_ALLOW_MINUS_ONE)

    def BeginEdit(self, row, col, grid):
        self.startValue = grid.GetTable().GetValue(row, col)
        self._tc.SetValue(str(self.startValue))
        self._tc.SetInsertionPointEnd()
        self._tc.SetFocus()
        self._tc.SetSelection(0, self._tc.GetLastPosition())

    def EndEdit(self, row, col, grid, oldVal=None):
        changed = False

        self.nextval = self.startValue
        val = self._tc.GetValue().lower().strip()
        val = val.replace(',','.') # allow , for decimal 
        if val != str(self.startValue):
            changed = True
            neg = False
            mult = False
            divide = False
            add = False
            if val.startswith('*'):
                mult = True
                val = val[1:]
            elif val.startswith('/'):
                divide = True
                val = val[1:]
            elif val.startswith('+'):
                add = True
                val = val[1:]
            if val.startswith('-'):
                neg = True
                val = val[1:]
            # allow old GSAS s20 and c20 etc for sind(20) and cosd(20)
            if val.startswith('s') and '(' not in val:
                val = 'sind('+val.strip('s')+')'
            elif val.startswith('c') and '(' not in val:
                val = 'cosd('+val.strip('c')+')'
            if neg:
                val = '-' + val
            val = G2fil.FormulaEval(val)
            if val is not None:
                if mult:
                    self.nextval *= val
                elif divide:
                    if val != 0: self.nextval /= val
                elif add:
                    self.nextval += val
                else:
                    self.nextval = val
            else:
                return None
            if oldVal is None: # this arg appears in 2.9+; before, we should go ahead & change the table
                grid.GetTable().SetValue(row, col, val) # update the table
            # otherwise self.ApplyEdit gets called

        self.startValue = ''
        self._tc.SetValue('')
        return changed
    
    def ApplyEdit(self, row, col, grid):
        """ Called only in wx >= 2.9
        Save the value of the control into the grid if EndEdit() returns as True
        """
        grid.GetTable().SetValue(row, col, self.nextval) # update the table

    def Reset(self):
        self._tc.SetValue(str(self.startValue))
        self._tc.SetInsertionPointEnd()

    def Clone(self,grid):
        return GridFractionEditor(grid)

    def StartingKey(self, evt):
        self.OnChar(evt)
        if evt.GetSkipped():
            self._tc.EmulateKeyPress(evt)

    def OnChar(self, evt):
        key = evt.GetKeyCode()
        if key < 32 or key >= 127: # outside printable ascii range; needed for backspace etc.
            evt.Skip()
        elif chr(key).lower() in '.+-*/0123456789cosind(),':
            evt.Skip()
        else:
            evt.StopPropagation()

#####  Get an output file or directory ################################################################################
def askSaveFile(G2frame,defnam,extension,longFormatName,parent=None):
    '''Ask the user to supply a file name; used for svn

    :param wx.Frame G2frame: The main GSAS-II window
    :param str defnam: a default file name
    :param str extension: the default file extension beginning with a '.'
    :param str longFormatName: a description of the type of file
    :param wx.Frame parent: the parent window for the dialog. Defaults
      to G2frame.

    :returns: a file name (str) or None if Cancel is pressed
    '''

    if not parent: parent = G2frame
    pth = GetExportPath(G2frame)
    #if GSASIIpath.GetConfigValue('debug'): print('debug: askSaveFile to '+pth)
    dlg = wx.FileDialog(
        parent, 'Input name for file to write', pth, defnam,
        longFormatName+' (*'+extension+')|*'+extension,
        wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
    dlg.CenterOnParent()
    try:
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            G2frame.LastExportDir = os.path.split(filename)[0]
            filename = os.path.splitext(filename)[0]+extension # make sure extension is correct
        else:
            filename = None
    finally:
        dlg.Destroy()
    return filename

def askSaveDirectory(G2frame):
    '''Ask the user to supply a directory name. Path name is used as the
    starting point for the next export path search. 

    :returns: a directory name (str) or None if Cancel is pressed
    '''
    pth = GetExportPath(G2frame)
    dlg = wx.DirDialog(G2frame,'Input directory where file(s) will be written',pth,wx.DD_DEFAULT_STYLE)
    dlg.CenterOnParent()
    try:
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            G2frame.LastExportDir = filename
        else:
            filename = None
    finally:
        dlg.Destroy()
    return filename

#####  Customized Notebook ################################################################################           
class GSNoteBook(wx.aui.AuiNotebook):
    '''Notebook used in various locations; implemented with wx.aui extension
    '''
    def __init__(self, parent, name='',size = None,style=wxaui_NB_TOPSCROLL):
        wx.aui.AuiNotebook.__init__(self, parent, style=style)
        if size: self.SetSize(size)
        self.parent = parent
        self.PageChangeHandler = None
        
    def PageChangeEvent(self,event):
        pass
                                                                  
    def Clear(self):
        GSNoteBook.DeleteAllPages(self)
        
    def FindPage(self,name):
        numPage = self.GetPageCount()
        for page in range(numPage):
            if self.GetPageText(page) == name:
                return page
        return None

    def ChangeSelection(self,page):
        # in wx.Notebook ChangeSelection is like SetSelection, but it
        # does not invoke the event related to pressing the tab button
        # I don't see a way to do that in aui.
        oldPage = self.GetSelection()
        self.SetSelection(page)
        return oldPage

    # def __getattribute__(self,name):
    #     '''This method provides a way to print out a message every time
    #     that a method in a class is called -- to see what all the calls
    #     might be, or where they might be coming from.
    #     Cute trick for debugging!
    #     '''
    #     attr = object.__getattribute__(self, name)
    #     if hasattr(attr, '__call__'):
    #         def newfunc(*args, **kwargs):
    #             print('GSauiNoteBook calling %s' %attr.__name__)
    #             result = attr(*args, **kwargs)
    #             return result
    #         return newfunc
    #     else:
    #         return attr
            
#### Help support routines ################################################################################
class MyHelp(wx.Menu):
    '''
    A class that creates the contents of a help menu.
    The menu will start with two entries:

    * 'Help on <helpType>': where helpType is a reference to an HTML page to
      be opened
    * About: opens an About dialog using OnHelpAbout. N.B. on the Mac this
      gets moved to the App menu to be consistent with Apple style.

    NOTE: for this to work properly with respect to system menus, the title
    for the menu must be &Help, or it will not be processed properly:

    ::

       menu.Append(menu=MyHelp(self,...),title="&Help")

    '''
    def __init__(self,frame,includeTree=False,morehelpitems=[]):
        wx.Menu.__init__(self,'')
        self.HelpById = {}
        self.frame = frame
        self.Append(wx.ID_ABOUT,'&About GSAS-II',
                        'Shows version and citation info')
        frame.Bind(wx.EVT_MENU, self.OnHelpAbout, id=wx.ID_ABOUT)
        if GSASIIpath.HowIsG2Installed():
            helpobj = self.Append(wx.ID_ANY,'&Check for updates\tCtrl+U',
                    'Updates to latest GSAS-II version')
            if os.access(GSASIIpath.path2GSAS2, os.W_OK):
                frame.Bind(wx.EVT_MENU, self.OnCheckUpdates, helpobj)
            else:
                helpobj.Enable(False)
            helpobj = self.Append(wx.ID_ANY,'&Regress to old GSAS-II version',
                    'Installs previous GSAS-II version')
            if os.access(GSASIIpath.path2GSAS2, os.W_OK):
                frame.Bind(wx.EVT_MENU, self.OnSelectVersion, helpobj)
            else:
                helpobj.Enable(False)
        if (GSASIIpath.HowIsG2Installed().startswith('git')
                and GSASIIpath.GetConfigValue('debug')): 
            helpobj = self.Append(wx.ID_ANY,'Switch to/from branch',
                    'Switch to/from a GSAS-II development branch')
            frame.Bind(wx.EVT_MENU, gitSelectBranch, helpobj)
        # provide special help topic names for extra items in help menu
        for lbl,indx in morehelpitems:
            helpobj = self.Append(wx.ID_ANY,lbl,'')
            frame.Bind(wx.EVT_MENU, self.OnHelpById, helpobj)
            self.HelpById[helpobj.GetId()] = indx
        # add help lookup(s) in gsasii.html
        self.AppendSeparator()
        if includeTree:
            helpobj = self.Append(wx.ID_ANY,'Help on GSAS-II',
                'Access web page with information on GSAS-II')
            frame.Bind(wx.EVT_MENU, self.OnHelpById, id=helpobj.GetId())
            self.HelpById[helpobj.GetId()] = 'Data tree'
        helpobj = self.Append(wx.ID_ANY,'Help on current data tree item\tF1',
                'Access web page on selected item in tree')
        frame.Bind(wx.EVT_MENU, self.OnHelpById, id=helpobj.GetId())
       
    def OnHelpById(self,event):
        '''Called when Help on... is pressed in a menu. Brings up a web page
        for documentation. Uses the helpKey value from the dataWindow window
        unless a special help key value has been defined for this menu id in
        self.HelpById

        Note that self should now (2frame) be child of the main window (G2frame)
        '''
        if hasattr(self.frame,'dataWindow'):  # Debug code: check this is called from menu in G2frame
            # should always be true in 2 Frame version
            dW = self.frame.dataWindow
        else:
            print('help error: not called from standard menu?')
            print (self)
            return            
        try:
            helpKey = dW.helpKey # look up help from helpKey in data window
            #if GSASIIpath.GetConfigValue('debug'): print 'DBG_dataWindow help: key=',helpKey
        except AttributeError:
            helpKey = ''
            if GSASIIpath.GetConfigValue('debug'): print('DBG_No helpKey for current dataWindow!')
        helpType = self.HelpById.get(event.GetId(),helpKey) # see if there is a special help topic
        #if GSASIIpath.GetConfigValue('debug'): print 'DBG_helpKey=',helpKey,'  helpType=',helpType
        if helpType == 'Tutorials':
            dlg = OpenTutorial(self.frame)
            dlg.ShowModal()
            dlg.Destroy()
            return
        else:
            ShowHelp(helpType,self.frame)

    def OnHelpAbout(self, event):
        "Display an 'About GSAS-II' box"
        try:
            import wx.adv as wxadv  # AboutBox moved here in Phoenix
        except:
            wxadv = wx
        info = wxadv.AboutDialogInfo()
        info.Name = 'GSAS-II'
        info.SetVersion(GSASIIpath.getG2VersionInfo())
        #info.Developers = ['Robert B. Von Dreele','Brian H. Toby']
        info.Copyright = ('(c) ' + time.strftime('%Y') +
''' Argonne National Laboratory
This product includes software developed
by the UChicago Argonne, LLC, as 
Operator of Argonne National Laboratory.''')
        info.Description = '''General Structure Analysis System-II (GSAS-II)
Robert B. Von Dreele and Brian H. Toby

Please cite as:
  B.H. Toby & R.B. Von Dreele, J. Appl. Cryst. 46, 544-549 (2013) 
For small angle use cite: 
  R.B. Von Dreele, J. Appl. Cryst. 47, 1748-9 (2014)
For DIFFaX use cite: 
  M.M.J. Treacy, J.M. Newsam & M.W. Deem, 
  Proc. Roy. Soc. Lond. A 433, 499-520 (1991)
'''
        info.WebSite = ("https://gsasii.github.io","GSAS-II home page")
        wxadv.AboutBox(info)

    def OnCheckUpdates(self,event):
        '''Check if the GSAS-II repository has an update for the current source files
        and perform that update if requested.
        '''
        if GSASIIpath.HowIsG2Installed().startswith('git'):
            gitCheckUpdates(self.frame)
        elif GSASIIpath.HowIsG2Installed().startswith('svn'):
            msg = '''
Note: GSAS-II is migrating to GitHub from the APS subversion server (where 
your current version of GSAS-II was installed from.) 

To use the new server, you must reinstall GSAS-II from https://bit.ly/G2download
(https://github.com/AdvancedPhotonSource/GSAS-II-buildtools/releases/tag/v1.0.1).
After September 2024, it will not be possible to update, unless you reinstall.

See web page GSASII.github.io for information on how to install.
'''
            res = ShowScrolledInfo(self.frame,msg,header='Please Note',
                                height=200,
                                buttonlist=[
           ('Open download site',
            lambda event: event.GetEventObject().GetParent().EndModal(wx.ID_OK)),
           ('Skip download for now',
            lambda event: event.GetEventObject().GetParent().EndModal(wx.ID_CANCEL))
                                    ])
            if res == wx.ID_OK:
                ShowWebPage(
                    'https://github.com/AdvancedPhotonSource/GSAS-II-buildtools/releases/tag/v1.0.1',
                    self.frame,
                    browser=True)
            svnCheckUpdates(self.frame)
        else:
            dlg = wx.MessageDialog(self.frame,
                                   'No VCS','Cannot update GSAS-II because it was not installed with a version control system or the VCS system could not be accessed.',
                                   wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return

    def OnSelectVersion(self,event):
        '''Allow the user to select a specific version of GSAS-II
        '''
        if GSASIIpath.HowIsG2Installed().startswith('git'):
            gitSelectVersion(self.frame)
        elif GSASIIpath.HowIsG2Installed().startswith('svn'):
            svnSelectVersion(self.frame)
        else:
            dlg = wx.MessageDialog(self.frame,
                                   'No VCS','Cannot update GSAS-II because it was not installed with a version control system or the VCS system could not be accessed.',
                                   wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return

################################################################################
class HelpButton(wx.Button):
    '''Create a help button that displays help information.
    The text can be displayed in a modal message window or it can be 
    a reference to a location in the gsasII.html (etc.) help web page, in which 
    case that page is opened in a web browser. 

    TODO: it might be nice if it were non-modal: e.g. it stays around until
    the parent is deleted or the user closes it, but this did not work for
    me. 

    :param parent: the panel/frame where the button will be placed
    :param str msg: the help text to be displayed. Indentation on 
       multiline help text is stripped (see :func:`StripIndents`). If wrap
       is set as non-zero, all new lines are 
    :param str helpIndex: location of the help information in the gsasII.html
      help file in the form of an anchor string. The URL will be 
      constructed from: location + gsasII.html + "#" + helpIndex
    :param int wrap: if specified, the text displayed is reformatted by
      wrapping it to fit in wrap pixels. Default is None which prevents 
      wrapping.
    '''
    def __init__(self,parent,msg='',helpIndex='',wrap=None):
        if sys.platform == "darwin": 
            wx.Button.__init__(self,parent,wx.ID_HELP)
        else:
            wx.Button.__init__(self,parent,wx.ID_ANY,'?',style=wx.BU_EXACTFIT)
        self.Bind(wx.EVT_BUTTON,self._onPress)
        if wrap:
            self.msg=StripIndents(msg,True)
        else:
            self.msg=StripIndents(msg)
        self.parent = parent
        self.helpIndex = helpIndex
        self.wrap = wrap
        self.msg = self.msg.replace(' & ',' && ')
    def _onClose(self,event):
        self.dlg.EndModal(wx.ID_CANCEL)
    def _onPress(self,event):
        'Respond to a button press by displaying the requested text'
        if self.helpIndex:
            ShowHelp(self.helpIndex,self.parent)
            return
        self.dlg = wx.Dialog(self.parent,wx.ID_ANY,'Help information', 
                        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        txt = wx.StaticText(self.dlg,wx.ID_ANY,self.msg)
        if self.wrap:
            txt.Wrap(self.wrap)
        mainSizer.Add(txt,1,wx.ALL|wx.EXPAND,10)
        txt.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(self.dlg, wx.ID_CLOSE) 
        btn.Bind(wx.EVT_BUTTON,self._onClose)
        btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        self.dlg.SetSizer(mainSizer)
        mainSizer.Fit(self.dlg)
        self.dlg.CenterOnParent()
        self.dlg.ShowModal()
        self.dlg.Destroy()
################################################################################
updateNoticeDict = {4919:True}  # example: {1234:True, 5000:False}
'''A dict with versions that should be noted. The value associated with the
tag is if all older projects should show the warning, or only the first 
to be opened. 
'''
def updateNotifier(G2frame,fileVersion):
    '''Posts an update notice when a a specially tagged GSAS-II version 
    is seen for the first time. Versions to be tagged are set in global
    updateNoticeDict; version info is found in file versioninfo.txt. 

    :param wx.Frame G2frame: GSAS-II main window
    :param int fileVersion: version of GSAS-II used to create the current
      .gpx file
    '''
    def tblLine(dlg,pixels=3):
        'place line in table'
        txtbox = wx.StaticText(dlg,wx.ID_ANY,'',size=(-1,pixels))
        txtbox.SetBackgroundColour(wx.Colour(0,0,0))
        tblSizer.Add(txtbox,0,wx.EXPAND)
        txtbox = wx.StaticText(dlg,wx.ID_ANY,'',size=(-1,pixels))
        txtbox.SetBackgroundColour(wx.Colour(0,0,0))
        tblSizer.Add(txtbox,0,wx.EXPAND)
    size = (700,500)
    rev = GSASIIpath.GetVersionNumber()
    try:
        int(rev)
    except:
        pass
    lastNotice = max(GSASIIpath.GetConfigValue('lastUpdateNotice',0),fileVersion)
    show = None               # first version number to show
    allProjects = False       # if True notice is shown for all projects, otherwise only once
    for key in updateNoticeDict:
        if updateNoticeDict[key]:
            if key >= fileVersion:
                if show is None:
                    show = fileVersion
                else:
                    show = min(show,fileVersion)
                allProjects = True
        else:
            if key >= lastNotice:
                if show is None:
                    show = lastNotice
                else:
                    show = min(show,lastNotice)
    if show is None: return

    filnam = os.path.join(GSASIIpath.path2GSAS2,'inputs','versioninfo.txt')
    if not os.path.exists(filnam):  # patch 3/2024 for svn dir organization
        filnam = os.path.join(GSASIIpath.path2GSAS2,'versioninfo.txt')    
    if not os.path.exists(filnam):
        print('Warning: file versioninfo.txt not found')
        return
    fp = open(filnam, 'r')
    vers = None
    noticeDict = {}
    for line in fp: 
        if line.strip().startswith('#'): continue
        if vers is not None:
            if len(line.strip()) == 0:
                vers = None
                continue
            else:
                noticeDict[vers] += ' '
                noticeDict[vers] += line.strip().replace('%%','\n')
        elif ':' in line:
            vers,tag = line.strip().split(':',1)
            try:
                vers = int(vers)
            except:
                continue
            noticeDict[vers] = tag.strip().replace('%%','\n')
    fp.close()
    for key in list(noticeDict.keys()):
        if key <= show: del noticeDict[key]
    if len(noticeDict) == 0: return

    dlg = wx.Dialog(G2frame,wx.ID_ANY,'Update notices',
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    sizer = wx.BoxSizer(wx.VERTICAL)
    txtbox = wx.StaticText(dlg,wx.ID_ANY,
            'Please read the notices below about major GSAS-II updates since'
            ' this project was last saved')
    sizer.Add(txtbox,0)
    txtbox.Wrap(size[0]-10)
    sizer.Add((10,10))
    panel = wxscroll.ScrolledPanel(dlg, wx.ID_ANY, size=(size[0]-20, size[1]))
    sizer.Add(panel,1,wx.EXPAND,1)
    tblSizer = wx.FlexGridSizer(0,2,5,10)
    tblLine(panel)
    txtbox = wx.StaticText(panel,wx.ID_ANY,'Version')
    tblSizer.Add(txtbox,0,wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL)
    txtbox = wx.StaticText(panel,wx.ID_ANY,'Notice')
    tblSizer.Add(txtbox,0,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    tblLine(panel)
    sizer.Add((10,10))
    btnsizer = wx.StdDialogButtonSizer()
    if allProjects: # will be shown again
        OKbtn = wx.Button(dlg, wx.ID_OK)
    else:
        OKbtn = wx.Button(dlg, wx.ID_OK,label='Record as seen')
        btn = wx.Button(dlg, wx.ID_CANCEL,label='Show again')
        btnsizer.AddButton(btn)
    OKbtn.SetDefault()
    btnsizer.AddButton(OKbtn)
    btnsizer.Realize()
    sizer.Add((-1,5))
    sizer.Add(btnsizer,0,wx.ALIGN_CENTER,50)
    for i,key in enumerate(sorted(noticeDict,reverse=True)):
        if i != 0: tblLine(panel,1)
        txtbox = wx.StaticText(panel,wx.ID_ANY,str(key))
        txtbox.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
        tblSizer.Add(txtbox,0,wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL)
        txtbox = wx.StaticText(panel,wx.ID_ANY,noticeDict[key])
        txtbox.Wrap(size[0]-110)
        txtbox.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
        tblSizer.Add(txtbox)
    tblLine(panel)
    panel.SetSizer(tblSizer)
    panel.SetAutoLayout(1)
    panel.SetupScrolling()
    dlg.SetSizer(sizer)
    sizer.Fit(dlg)
    dlg.CenterOnParent()
    if dlg.ShowModal() == wx.ID_OK:
        if not allProjects:
            vars = GetConfigValsDocs()
            vars['lastUpdateNotice'][1] = rev
            GSASIIpath.SetConfigValue(vars)
            SaveConfigVars(vars)
        dlg.Destroy()
        return
    dlg.Destroy()
                  
################################################################################
class MyHtmlPanel(wx.Panel):
    '''Defines a panel to display HTML help information, as an alternative to
    displaying help information in a web browser.
    '''
    def __init__(self, frame, newId):
        self.frame = frame
        wx.Panel.__init__(self, frame, newId)
        sizer = wx.BoxSizer(wx.VERTICAL)
        back = wx.Button(self, -1, "Back")
        back.Bind(wx.EVT_BUTTON, self.OnBack)
        self.htmlwin = G2HtmlWindow(self, newId, size=(750,450))
        sizer.Add(self.htmlwin, 1,wx.EXPAND)
        sizer.Add(back, 0, wx.ALIGN_LEFT, 0)
        self.SetSizer(sizer)
        sizer.Fit(frame)        
        self.Bind(wx.EVT_SIZE,self.OnHelpSize)
    def OnHelpSize(self,event):         #does the job but weirdly!!
        anchor = self.htmlwin.GetOpenedAnchor()
        if anchor:            
            self.htmlwin.ScrollToAnchor(anchor)
            wx.CallAfter(self.htmlwin.ScrollToAnchor,anchor)
            if event: event.Skip()
    def OnBack(self, event):
        self.htmlwin.HistoryBack()
    def LoadFile(self,file):
        pos = file.rfind('#')
        if pos != -1:
            helpfile = file[:pos]
            helpanchor = file[pos+1:]
        else:
            helpfile = file
            helpanchor = None
        self.htmlwin.LoadPage(helpfile)
        if helpanchor is not None:
            self.htmlwin.ScrollToAnchor(helpanchor)
            xs,ys = self.htmlwin.GetViewStart()
            self.htmlwin.Scroll(xs,ys-1)
################################################################################
class G2HtmlWindow(wx.html.HtmlWindow):
    '''Displays help information in a primitive HTML browser type window
    '''
    def __init__(self, parent, *args, **kwargs):
        self.parent = parent
        wx.html.HtmlWindow.__init__(self, parent, *args, **kwargs)
    def LoadPage(self, *args, **kwargs):
        wx.html.HtmlWindow.LoadPage(self, *args, **kwargs)
        self.TitlePage()
    def OnLinkClicked(self, *args, **kwargs):
        wx.html.HtmlWindow.OnLinkClicked(self, *args, **kwargs)
        xs,ys = self.GetViewStart()
        self.Scroll(xs,ys-1)
        self.TitlePage()
    def HistoryBack(self, *args, **kwargs):
        wx.html.HtmlWindow.HistoryBack(self, *args, **kwargs)
        self.TitlePage()
    def TitlePage(self):
        self.parent.frame.SetTitle(self.GetOpenedPage() + ' -- ' + 
            self.GetOpenedPageTitle())

################################################################################
def StripIndents(msg,singleLine=False):
    '''Strip unintended indentation from multiline strings. 
    When singleLine is True, all newline are removed, but inserting "%%"
    into the string will cause a blank line to be inserted at that point
    and %t% will generate a new line and tab (to indent a line)

    :param str msg: a string containing one or more lines of text. 
      spaces or tabs following a newline are removed.
    :param bool singleLine: removes all newlines from the msg so that 
      the text may be wrapped. 
    :returns: the string but reformatted
    '''
    msg1 = msg.replace('\n ','\n')
    while msg != msg1:
        msg = msg1
        msg1 = msg.replace('\n ','\n')
        msg1 = msg1.replace('  ',' ')
    msg1 = msg.replace('\n\t','\n')
    while msg != msg1:
        msg = msg1
        msg1 = msg.replace('\n\t','\n')
    if singleLine:
        msg = msg.replace('\n',' ')
        msg1 = msg.replace('  ',' ')
        while msg != msg1:
            msg = msg1
            msg1 = msg1.replace('  ',' ')
        msg = msg.replace('%%','\n\n')
        msg = msg.replace('%t%','\n\t')
    return msg

def StripUnicode(string,subs='.'):
    '''Strip non-ASCII characters from strings
    
    :param str string: string to strip Unicode characters from
    :param str subs: character(s) to place into string in place of each
      Unicode character. Defaults to '.'

    :returns: a new string with only ASCII characters
    '''
    s = ''
    for c in string:
        if ord(c) < 128:
            s += c
        else:
            s += subs
    return s.encode('ascii','replace')

def getTextSize(txt):
    'Get the size of the text string txt in points, returns (x,y)'
    dc = wx.ScreenDC()
    return tuple(dc.GetTextExtent(txt))
    
# wx classes for reading various types of data files ######################################################################
def BlockSelector(ChoiceList, ParentFrame=None,title='Select a block',
    size=None, header='Block Selector',useCancel=True):
    ''' Provide a wx dialog to select a single block where one must 
    be selected. Used for selecting for banks for instrument 
    parameters if the file contains more than one set.
    '''
    if useCancel:
        dlg = wx.SingleChoiceDialog(
            ParentFrame,title, header,ChoiceList)
    else:
        dlg = wx.SingleChoiceDialog(
            ParentFrame,title, header,ChoiceList,
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.OK|wx.CENTRE)
    if size: dlg.SetMinSize(size)
    dlg.CenterOnParent()
    dlg.SendSizeEvent()
    if dlg.ShowModal() == wx.ID_OK:
        sel = dlg.GetSelection()
        return sel
    else:
        return None
    dlg.Destroy()

def MultipleBlockSelector(ChoiceList, ParentFrame=None,
    title='Select a block',size=None, header='Block Selector'):
    '''Provide a wx dialog to select a block of data if the
    file contains more than one set of data and one must be
    selected. Used in :mod:`G2pwd_CIF` only.

    :returns: a list of the selected blocks
    '''
    dlg = wx.MultiChoiceDialog(ParentFrame,title, header,ChoiceList+['Select all'],
        wx.CHOICEDLG_STYLE)
    if size: dlg.SetMinSize(size)
    dlg.CenterOnScreen()
    dlg.SendSizeEvent()
    if dlg.ShowModal() == wx.ID_OK:
        sel = dlg.GetSelections()
    else:
        return []
    dlg.Destroy()
    selected = []
    if len(ChoiceList) in sel:
        return range(len(ChoiceList))
    else:
        return sel
    return selected

class MultipleChoicesDialog(wx.Dialog):
    '''A dialog that offers a series of choices, each with a
    title and a wx.Choice widget. Intended to be used Modally. 
    typical input:

        *  choicelist=[ ('a','b','c'), ('test1','test2'),('no choice',)]
        *  headinglist = [ 'select a, b or c', 'select 1 of 2', 'No option here']
        
    selections are placed in self.chosen when OK is pressed

    Also see GSASIIctrlGUI
    '''
    def __init__(self,choicelist,headinglist,
                 head='Select options',
                 title='Please select from options below',
                 parent=None):
        self.chosen = []
        wx.Dialog.__init__(
            self,parent,wx.ID_ANY,head, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((10,10),1)
        topLabl = wx.StaticText(panel,wx.ID_ANY,title)
        mainSizer.Add(topLabl,0,wx.ALIGN_CENTER_VERTICAL|wx.CENTER,10)
        self.ChItems = []
        for choice,lbl in zip(choicelist,headinglist):
            mainSizer.Add((10,10),1)
            self.chosen.append(0)
            topLabl = wx.StaticText(panel,wx.ID_ANY,' '+lbl)
            mainSizer.Add(topLabl,0,wx.ALIGN_LEFT,10)
            self.ChItems.append(wx.Choice(self, wx.ID_ANY, (100, 50), choices = choice))
            mainSizer.Add(self.ChItems[-1],0,wx.ALIGN_CENTER,10)

        OkBtn = wx.Button(panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        cancelBtn = wx.Button(panel,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        panel.SetSizer(mainSizer)
        panel.Fit()
        self.Fit()
        
    def OnOk(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        # save the results from the choice widgets
        self.chosen = []
        for w in self.ChItems:
            self.chosen.append(w.GetSelection())
        self.EndModal(wx.ID_OK)              
            
    def OnCancel(self,event):
        parent = self.GetParent()
        if parent is not None: parent.Raise()
        self.chosen = []
        self.EndModal(wx.ID_CANCEL)              

def MultipleChoicesSelector(choicelist, headinglist, ParentFrame=None, **kwargs):
    '''A modal dialog that offers a series of choices, each with a title and a wx.Choice
    widget. Used in :mod:`G2pwd_CIF` only.

    Typical input:
    
       * choicelist=[ ('a','b','c'), ('test1','test2'),('no choice',)]
       
       * headinglist = [ 'select a, b or c', 'select 1 of 2', 'No option here']
       
    optional keyword parameters are: head (window title) and title
    returns a list of selected indicies for each choice (or None)
    '''
    result = None
    dlg = MultipleChoicesDialog(choicelist,headinglist,
        parent=ParentFrame, **kwargs)          
    dlg.CenterOnParent()
    if dlg.ShowModal() == wx.ID_OK:
        result = dlg.chosen
    dlg.Destroy()
    return result

def PhaseSelector(ChoiceList, ParentFrame=None,
    title='Select a phase', size=None,header='Phase Selector'):
    ''' Provide a wx dialog to select a phase, used in importers if a file 
    contains more than one phase
    '''
    return BlockSelector(ChoiceList,ParentFrame,title,
        size,header)

def XformMatrix(panel,Trans,Uvec,Vvec,OnLeave=None,OnLeaveArgs={}):
    '''Display a transformation matrix and two vectors'''
    Trmat = wx.FlexGridSizer(4,6,0,0)
    Trmat.Add((10,0),0)
    Trmat.Add(wx.StaticText(panel,label='      M'),wx.ALIGN_CENTER)
    Trmat.Add((10,0),0)
    Trmat.Add((10,0),0)
    Trmat.Add(wx.StaticText(panel,label='      U'),wx.ALIGN_CENTER)
    Trmat.Add(wx.StaticText(panel,label='      V'),wx.ALIGN_CENTER)
        
    for iy,line in enumerate(Trans):
        for ix,val in enumerate(line):
            item = ValidatedTxtCtrl(panel,Trans[iy],ix,nDig=(10,3),size=(65,25),
                                    OnLeave=OnLeave,OnLeaveArgs=OnLeaveArgs)
            Trmat.Add(item)
        Trmat.Add((25,0),0)
        vec = ValidatedTxtCtrl(panel,Uvec,iy,nDig=(10,3),size=(65,25),
                                    OnLeave=OnLeave,OnLeaveArgs=OnLeaveArgs)
        Trmat.Add(vec)
        vec = ValidatedTxtCtrl(panel,Vvec,iy,nDig=(10,3),size=(65,25),
                                    OnLeave=OnLeave,OnLeaveArgs=OnLeaveArgs)
        Trmat.Add(vec)
    return Trmat

def showUniqueCell(frame,cellSizer,row,cell,SGData=None,
                       editAllowed=False,OnCellChange=None):
    '''function to put cell values into a GridBagSizer. 
    First column (#0) is reserved for labels etc. 
    if editAllowed is True, values are placed in a wx.TextCtrl and if needed
    two rows are used in the table.
    '''
    cellGUIlist = [
        [['m3','m3m'],[" Unit cell: a = "],["{:.5f}"],[0]],
        [['3R','3mR'],[" a = ",u" \u03B1 = "],["{:.5f}","{:.3f}"],[0,3]],
        [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],
             [" a = "," c = ",u" \u03B3 = "],
             ["{:.5f}","{:.5f}","{:.3f}"],[0,2,-5]],
        [['mmm'],[" a = "," b = "," c = "],["{:.5f}","{:.5f}","{:.5f}"],
            [0,1,2]],
        [['2/m'+'a'],[" a = "," b = "," c = ",u" \u03B1 = "],
            ["{:.5f}","{:.5f}","{:.5f}","{:.3f}"],[0,1,2,3]],
        [['2/m'+'b'],[" a = "," b = "," c = ",u" \u03B2 = "],
            ["{:.5f}","{:.5f}","{:.5f}","{:.3f}"],[0,1,2,4]],
        [['2/m'+'c'],[" a = "," b = "," c = ",u" \u03B3 = "],
            ["{:.5f}","{:.5f}","{:.5f}","{:.3f}"],[0,1,2,5]],
        [['-1'],[" a = "," b = "," c = ",u" \u03B1 = ",u" \u03B2 = ",u" \u03B3 = "],
             ["{:.5f}","{:.5f}","{:.5f}","{:.3f}","{:.3f}","{:.3f}"],[0,1,2,3,4,5]]
        ]
    cellList = []
    if SGData is None:
        laue = '-1'
    else:
        laue = SGData['SGLaue']
        if laue == '2/m': laue += SGData['SGUniq']
    for cellGUI in cellGUIlist:
        if laue in cellGUI[0]:
            useGUI = cellGUI
            break
    for txt,fmt,indx in zip(*useGUI[1:]):
        col = 1+2*abs(indx)
        cellrow = row
        if editAllowed and abs(indx) > 2:
            cellrow = row + 1
            col = 1+2*(abs(indx)-3)
        cellSizer.Add(wx.StaticText(frame,label=txt),(cellrow,col))
        if editAllowed and indx >= 0:
            Fmt = (10,5)
            if '.3' in fmt: Fmt = (10,3) 
            cellVal = ValidatedTxtCtrl(frame,cell,indx,
                    xmin=0.1,xmax=500.,nDig=Fmt,OnLeave=OnCellChange)
            cellSizer.Add(cellVal,(cellrow,col+1))
            cellList.append(cellVal.GetId())
        else:
            cellSizer.Add(wx.StaticText(frame,label=fmt.format(cell[abs(indx)])),(cellrow,col+1))
    #volume
    volCol = 13
    if editAllowed: 
        volCol = 8
    cellSizer.Add(wx.StaticText(frame,label=' Vol = '),(row,volCol))
    if editAllowed:
        volVal = wx.TextCtrl(frame,value=('{:.2f}'.format(cell[6])),style=wx.TE_READONLY)
        volVal.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE))
        cellSizer.Add(volVal,(row,volCol+1))
    else:
        cellSizer.Add(wx.StaticText(frame,label='{:.2f}'.format(cell[6])),(row,volCol+1))
    return cellrow,cellList


################################################################################
# configuration routines (for editing config.py)
def SaveGPXdirectory(path,write=True):
    if GSASIIpath.GetConfigValue('Starting_directory') == path: return
    vars = GetConfigValsDocs()
    try:
        vars['Starting_directory'][1] = path
        if GSASIIpath.GetConfigValue('debug'): print('DBG_Saving GPX path: '+path)
        if write: SaveConfigVars(vars)
    except KeyError:
        pass

def SaveImportDirectory(path):
    if GSASIIpath.GetConfigValue('Import_directory') == path: return
    vars = GetConfigValsDocs()
    try:
        vars['Import_directory'][1] = path
        if GSASIIpath.GetConfigValue('debug'): print('DBG_Saving Import path: '+path)
        SaveConfigVars(vars)
    except KeyError:
        pass

def GetConfigValsDocs():
    '''Reads the module referenced in fname (often <module>.__file__) and
    return a dict with names of global variables as keys.
    For each global variable, the value contains four items:

    :returns: a dict where keys are names defined in module config_example.py
      where the value is a list of four items, as follows:

         * item 0: the default value
         * item 1: the current value
         * item 2: the initial value (starts same as item 1)
         * item 3: the "docstring" that follows variable definition

    '''
    import config_example
    import ast
    fname = os.path.splitext(config_example.__file__)[0]+'.py' # convert .pyc to .py
    if '3' in platform.python_version_tuple()[0]: 
        with open(fname, 'r',encoding='utf-8') as f:
            fstr = f.read()
    else:
        with open(fname, 'r') as f:
            fstr = f.read()
    fstr = fstr.replace('\r\n', '\n').replace('\r', '\n')
    if not fstr.endswith('\n'):
        fstr += '\n'
    tree = ast.parse(fstr)
    d = {}
    key = None
    for node in ast.walk(tree):
        if isinstance(node,ast.Assign):
            key = node.targets[0].id
            d[key] = [config_example.__dict__.get(key),
                      GSASIIpath.configDict.get(key),
                      GSASIIpath.configDict.get(key),'']
        elif isinstance(node,ast.Expr) and key:
            d[key][3] = node.value.s.strip()
        else:
            key = None
    return d

inhibitSave = False
def SaveConfigVars(vars,parent=None):
    '''Write the current config variable values to config.py

    :params dict vars: a dictionary of variable settings and meanings as
      created in :func:`GetConfigValsDocs`.
    :param parent: wx.Frame object or None (default) for parent
      of error message if no file can be written.
    :returns: True if unable to write the file, None otherwise
    '''
    if inhibitSave:
        if GSASIIpath.GetConfigValue('debug'):
            print('inhibitSave prevents saving configuration')
        return

    # try to write to where an old config file is located
    try:
        import config
        savefile = config.__file__
    except ImportError: # no config.py file yet
        savefile = os.path.join(GSASIIpath.path2GSAS2,'config.py')
    except Exception: # import failed
        # find the bad file, save it in a new name and prepare to overwrite it
        for p in sys.path:
            savefile = os.path.join(p,'config.py')
            if os.path.exists(savefile):
                import distutils.file_util as dfu
                keepfile = os.path.join(p,'config.py_corrupt')
                print('Current config file contains an error:',savefile)
                print('saving that file as',keepfile)
                dfu.copy_file(savefile,keepfile)
                print('preparing to overwrite...')
                break
        else:
            print('unexpected error importing config.py')
            savefile = os.path.join(GSASIIpath.path2GSAS2,'config.py')
        
    # try to open file for write
    try:
        savefile = os.path.splitext(savefile)[0]+'.py' # convert .pyc to .py
        fp = open(savefile,'w')
    except IOError:  # can't write there, write in local mods directory
        # create a local mods directory, if needed
        g2local = os.path.expanduser('~/.G2local/')
        if not os.path.exists(g2local):
            try:
                print(u'Creating directory '+g2local)
                os.mkdir(g2local)
            except:
                if parent:
                    G2MessageBox(parent,
                                     f'Error trying to create directory {g2local}. Unable to save')
                else:
                    print(u'Error trying to create directory '+g2local)
                return True
            sys.path.insert(0,os.path.expanduser('~/.G2local/'))
        savefile = os.path.join(os.path.expanduser('~/.G2local/'),'config.py')
        try:
            fp = open(savefile,'w')
        except IOError:
            if parent:
                G2MessageBox(parent,'Error trying to write configuration to '+savefile,
                    'Unable to save')
            else:
                print('Error trying to write configuration to '+savefile)
            return True
    import datetime
    fp.write("# -*- coding: utf-8 -*-\n'''\n")
    fp.write("*config.py: Configuration options*\n----------------------------------\n")
    fp.write("This file created in SelectConfigSetting on {:%d %m %Y %H:%M}\n".
             format(datetime.datetime.now()))
    fp.write("'''\n\n")
    fp.write("import os.path\n")
    fp.write("import GSASIIpath\n\n")
    for var in sorted(vars.keys(),key=lambda s: s.lower()):
        if vars[var][1] is None: continue
        if vars[var][1] == '': continue
        if vars[var][0] == vars[var][1]: continue
        try:
            float(vars[var][1]) # test for number
            fp.write(var + ' = ' + str(vars[var][1])+'\n')
        except:
            if type(vars[var][1]) is list:
                fp.write(var + ' = [\n')
                for varstr in vars[var][1]:
                    if '\\' in varstr:
                        fp.write('\t  os.path.normpath("' + varstr.replace('\\','/') +'"),\n')
                    else:
                        fp.write('\t  "' + str(varstr)+'",\n')
                fp.write('   ]\n')
            elif type(vars[var][1]) is tuple:
                fp.write(var + ' = ' + str(vars[var][1])+'\n')
            else:
                try:
                    eval(vars[var][1]) # test for an expression
                    fp.write(var + ' = ' + str(vars[var][1])+'\n')
                except: # must be a string
                    varstr = vars[var][1]
                    if '\\' in varstr:
                        fp.write(var + ' = os.path.normpath("' + varstr.replace('\\','/') +'")\n')
                    else:
                        fp.write(var + ' = "' + str(varstr)+'"\n')
        if vars[var][3]:
            fp.write("'''" + str(vars[var][3]) + "\n'''\n\n")
    fp.close()
    print('wrote file '+savefile)
    
class SelectConfigSetting(wx.Dialog):
    '''Dialog to select configuration variables and set associated values.
    '''
    def __init__(self,parent):
        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, wx.ID_ANY, 'Set Config Variable', style=style)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.vars = GetConfigValsDocs()
        self.G2frame = parent
        self.restart = False

        label = wx.StaticText(
            self,  wx.ID_ANY,
            'Select a GSAS-II configuration variable to change'
            )
        self.sizer.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        self.choice = {}
        choices = sorted([k for k in self.vars if not k.startswith('enum_')],
                         key=lambda s: s.lower())
        btn = G2ChoiceButton(self, choices,
            strLoc=self.choice,strKey=0,onChoice=self.OnSelection)
        btn.SetLabel("")
        self.sizer.Add(btn)

        self.varsizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.varsizer,1,wx.ALL|wx.EXPAND,1)
        
        self.doclbl = wx.StaticBox(self, wx.ID_ANY, "")
        self.doclblsizr = wx.StaticBoxSizer(self.doclbl)
        self.docinfo = wx.StaticText(self,  wx.ID_ANY, "")
        self.doclblsizr.Add(self.docinfo, 0, wx.ALIGN_LEFT|wx.ALL, 5)
        self.sizer.Add(self.doclblsizr, 0, wx.EXPAND|wx.ALL, 5)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.saveBtn = wx.Button(self,-1,"Save current settings")
        btnsizer.Add(self.saveBtn, 0, wx.ALL, 2) 
        self.saveBtn.Bind(wx.EVT_BUTTON, self.OnSave)
        self.saveBtn.Enable(False)
        
        btn = wx.Button(self,wx.ID_CANCEL)
        btnsizer.Add(btn, 0, wx.ALL, 2) 
        self.sizer.Add(btnsizer, 0, wx.ALIGN_CENTRE|wx.ALL, 5) 
                
        self.SetSizer(self.sizer)
        self.sizer.Fit(self)
        self.CenterOnParent()
        
    def OnChange(self,event=None):
        ''' Check if anything been changed. Turn the save button on/off.
        '''
        for var in self.vars:
            if self.vars[var][0] is None and self.vars[var][1] is not None:
                # make blank strings into None, if that is the default
                if self.vars[var][1].strip() == '': self.vars[var][1] = None
            if self.vars[var][1] != self.vars[var][2]:
                #print 'changed',var,self.vars[var][:3]
                self.saveBtn.Enable(True)
                if 'restart' in self.vars[var][3].lower():
                    self.restart  = True
                break
        else:
            self.saveBtn.Enable(False)
        try:
            self.resetBtn.Enable(True)
        except:
            pass
        
    def OnApplyChanges(self,event=None):
        'Set config variables to match the current settings'
        GSASIIpath.SetConfigValue(self.vars)
        self.EndModal(wx.ID_OK)
        global inhibitSave
        if event is not None: inhibitSave = True
        import GSASIImpsubs as G2mp
        G2mp.ResetMP()
        
    def OnSave(self,event):
        '''Write the config variables to config.py and then set them
        as the current settings
        '''
        global inhibitSave
        inhibitSave = False
        if not SaveConfigVars(self.vars,parent=self):
            self.OnApplyChanges() # force a reload of the config settings
        else:
            self.EndModal(wx.ID_OK)

    def OnBoolSelect(self,event):
        'Respond to a change in a True/False variable'
        rb = event.GetEventObject()
        var = self.choice[0]
        self.vars[var][1] = (rb.GetSelection() == 0)
        self.OnChange()
        wx.CallAfter(self.OnSelection)
        
    def onSelDir(self,event):
        'Select a directory from a menu'
        dlg = wx.DirDialog(self, "Choose a directory:",style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            var = self.choice[0]
            self.vars[var][1] = dlg.GetPath()
            self.strEd.SetValue(self.vars[var][1])
            self.OnChange()
        dlg.Destroy()

    def onSelExec(self,event):
        'Select an executable file from a menu'
        var = self.choice[0]
        is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)
        defD = defF = ''
        if self.vars[var][1] is not None and os.path.exists(self.vars[var][1]):
            defD,defF=os.path.split(self.vars[var][1])
        repeat = True
        while repeat:
            repeat = False
            if sys.platform == "win32":
                dlg = wx.FileDialog(self, "Choose a .exe file:",
                                defaultDir=defD,defaultFile=defF,
                                style=wx.FD_DEFAULT_STYLE|wx.FD_FILE_MUST_EXIST,
                                wildcard="Executable files|*.exe")
            else:
                dlg = wx.FileDialog(self, "Choose an executable image:",
                                defaultDir=defD,defaultFile=defF,
                                style=wx.FD_DEFAULT_STYLE|wx.FD_FILE_MUST_EXIST)
            if dlg.ShowModal() == wx.ID_OK:
                val = dlg.GetPath()
                if os.path.exists(val) and is_exe(val):
                    self.vars[var][1] = val
                    self.strEd.SetValue(self.vars[var][1])
                    self.OnChange()
                else:
                    dlg.Destroy()
                    G2MessageBox(self,'File not found or not executable',
                                     'Invalid file')
                    repeat = True
                    continue
        dlg.Destroy()

        
    def OnSelection(self):
        'show a selected variable and allow it to be changed'
        def OnNewColorBar(event):
            self.vars['Contour_color'][1] = self.colSel.GetValue()
            self.OnChange(event)

        if 'phoenix' in wx.version():
            self.varsizer.Clear(True)
        else:
            self.varsizer.DeleteWindows()
        var = self.choice[0]
        showdef = True
        if var not in self.vars:
            raise Exception("How did this happen?")
        if 'enum_'+var in self.vars:
            choices = self.vars['enum_'+var][0]
            self.colSel = EnumSelector(self,self.vars[var],1,choices,
                                               OnChange=self.OnChange)       
            self.varsizer.Add(self.colSel, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        elif type(self.vars[var][0]) is int:
            ed = ValidatedTxtCtrl(self,self.vars[var],1,typeHint=int,OKcontrol=self.OnChange)
            self.varsizer.Add(ed, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        elif type(self.vars[var][0]) is float:
            ed = ValidatedTxtCtrl(self,self.vars[var],1,typeHint=float,OKcontrol=self.OnChange)
            self.varsizer.Add(ed, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        elif type(self.vars[var][0]) is bool:
            showdef = False
            lbl = "value for "+var
            ch = []
            for i,v in enumerate((True,False)):
                s = str(v)
                if v == self.vars[var][0]:
                    defopt = i
                    s += ' (default)'
                ch += [s]
            rb = wx.RadioBox(self, wx.ID_ANY, lbl, wx.DefaultPosition, wx.DefaultSize,
                ch, 1, wx.RA_SPECIFY_COLS)
            # set initial value
            if self.vars[var][1] is None:
                rb.SetSelection(defopt)
            elif self.vars[var][1]:
                rb.SetSelection(0)
            else:
                rb.SetSelection(1)
            rb.Bind(wx.EVT_RADIOBOX,self.OnBoolSelect)
            self.varsizer.Add(rb, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        else:
            if var.endswith('_directory') or var.endswith('_location'):
                btn = wx.Button(self,wx.ID_ANY,'Select from dialog...')
                btn.Bind(wx.EVT_BUTTON,self.onSelDir)
                sz = (400,-1)
            elif var.endswith('_exec'):
                btn = wx.Button(self,wx.ID_ANY,'Select from dialog...')
                btn.Bind(wx.EVT_BUTTON,self.onSelExec)
                sz = (400,-1)
            else:
                btn = None
                sz = (250,-1)
            if var == 'Contour_color':
                if self.vars[var][1] is None:
                    self.vars[var][1] = 'Paired'
                colorList = sorted([m for m in mpl.cm.datad.keys()]+['GSPaired','GSPaired_r',],key=lambda s: s.lower())   #if not m.endswith("_r")
                self.colSel = wx.ComboBox(self,value=self.vars[var][1],choices=colorList,
                    style=wx.CB_READONLY|wx.CB_DROPDOWN)
                self.colSel.Bind(wx.EVT_COMBOBOX, OnNewColorBar)
                self.varsizer.Add(self.colSel, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
            elif var == 'Image_calibrant':
                import ImageCalibrants as calFile
                calList = sorted([m for m in calFile.Calibrants.keys()],
                                     key=lambda s: s.lower())
                self.colSel = EnumSelector(self,self.vars[var],1,calList,
                                               OnChange=self.OnChange)       
                self.varsizer.Add(self.colSel, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
            else:
                self.strEd = ValidatedTxtCtrl(self,self.vars[var],1,typeHint=str,
                    OKcontrol=self.OnChange,size=sz)
                if self.vars[var][1] is not None:
                    self.strEd.SetValue(self.vars[var][1])
                if btn:
                    self.varsizer.Add(self.strEd, 0, wx.ALL|wx.EXPAND, 5)
                    self.varsizer.Add(btn, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
                else:
                    self.varsizer.Add(self.strEd, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        # button for reset to default value
        lbl = "Reset to Default"
        if showdef: # spell out default when needed
            lbl += ' (='+str(self.vars[var][0])+')'
            #label = wx.StaticText(self,  wx.ID_ANY, 'Default value = '+str(self.vars[var][0]))
            #self.varsizer.Add(label, 0, wx.ALIGN_LEFT|wx.ALL, 5)
        self.resetBtn = wx.Button(self,-1,lbl)
        self.resetBtn.Bind(wx.EVT_BUTTON, self.OnClear)
        if self.vars[var][1] is not None and self.vars[var][1] != '': # show current value, if one
            #label = wx.StaticText(self,  wx.ID_ANY, 'Current value = '+str(self.vars[var][1]))
            #self.varsizer.Add(label, 0, wx.ALIGN_LEFT|wx.ALL, 5)
            self.resetBtn.Enable(True)
        else:
            self.resetBtn.Enable(False)
        self.varsizer.Add(self.resetBtn, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        # show meaning, if defined
        self.doclbl.SetLabel("Description of "+str(var))
        vartxt = self.vars[var][3]
        vartxt = StripIndents(vartxt.replace(' & ',' && '),True)
        if vartxt:
            self.docinfo.SetLabel(vartxt)
        else:
            self.docinfo.SetLabel("(not documented)")
        self.docinfo.Wrap(500)
        self.sizer.Fit(self)
        self.CenterOnParent()
        wx.CallAfter(self.SendSizeEvent)

    def OnClear(self, event):
        var = self.choice[0]
        self.vars[var][1] = self.vars[var][0]
        self.OnChange()
        wx.CallAfter(self.OnSelection)

################################################################################
class RefinementProgress(wx.ProgressDialog):
    '''Defines a wrapper to place around wx.ProgressDialog to be used for 
    showing refinement progress. At some point a better progress window should be
    created that keeps useful info on the screen such as some starting and  
    current fit metrics, but for now all this adds is window defaults 
    and a wx.Yield call during progress update calls.
    '''
    def __init__(self, title='Residual', message='All data Rw =',
                     maximum=101, parent=None, trialMode=False, seqLen=0,
                     style=None):
        if style is None:
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT
        super(self.__class__,self).__init__(title, message, int(maximum), parent, style)
        Size = self.GetSize()
        if 50 < Size[0] < 500: # sanity check on size, since this fails w/Win & wx3.0
            self.SetSize((int(Size[0]*1.2),Size[1])) # increase size a bit along x
        self.CenterOnParent()
        self.Show()
        
    def Update(self,value, newmsg=""):
        wx.GetApp().Yield()
        #print('wx Yield called')
        #print('Update:',value,newmsg)
        return super(self.__class__,self).Update(int(value), newmsg)
        
################################################################################
fmtRw = lambda value: '{:.2f}'.format(float(value))
class G2RefinementProgress(wx.Dialog):
    '''Defines an replacement for wx.ProgressDialog to be used for 
    showing refinement progress.

    :param str title: string to place on border of window (default is 
            'Refinement progress').
    :param str message: initial string to place on top line of window.
    :param int maximum: maximum value for progress gauge bar on bottom 
       of window. 
    :param wx.Frame parent: parent window for creation of this dialog
    :param bool trialMode: Set to True for Levenberg-Marquardt fitting
      where Rw may be computed several times for a single cycle.
      Call :meth:`AdvanceCycle` when trialMode is True to indicate that a cycle
      has been completed. Default is False.
    :param int seqLen: Number of histograms in sequential fit. A value of 
      zero (default) means that the fit is not a sequential fit.
    :param int seqShow: Number of histograms to shown in a sequential fit (default 3)
    :param int style: optional parameters that determine how the dialog is
      is displayed.  
    '''
    def __init__(self, title='Refinement progress', message='All data Rw =',
        maximum=101, parent=None, trialMode=False,seqLen=0, seqShow=3,style=None):

        self.trialRw = trialMode # used for Levenberg-Marquardt fitting
        self.SeqLen = seqLen
        self.seqShow = seqShow
        if self.SeqLen:
            self.maxCycle = self.SeqLen        
        self.SeqCount = -1
        self.rows = 4
        if self.trialRw: self.rows = 5
        if style is None: style = wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER

        super(self.__class__,self).__init__(parent, wx.ID_ANY, title,style=style, size=(-1,-1))
        self.Bind(wx.EVT_CLOSE, self._onClose)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.startTime = time.time()
        self.abortStatus = False
        self.SetSizer(mainSizer)
        mainSizer.Add((-1,3))
        self.msgLine1 = wx.StaticText(self,wx.ID_ANY,message)
        mainSizer.Add(self.msgLine1)
        mainSizer.Add((-1,3))
        self.msgLine2 = wx.StaticText(self,wx.ID_ANY,'')
        mainSizer.Add(self.msgLine2)

        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        vSizer = wx.BoxSizer(wx.VERTICAL)
        lblSizer = self._makeLabeledTable()   # make the Table for Rfactors in the RwPanel
        vSizer.Add((-1,-1),1,wx.EXPAND,1)
        vSizer.Add(lblSizer,0,wx.EXPAND)
        vSizer.Add((-1,-1),1,wx.EXPAND,1)
        hSizer.Add(vSizer,1,wx.EXPAND,1)
        pltPanel = wx.Panel(self,size=(-1,-1))
        self.figure = mplfig.Figure(dpi=100,figsize=(3,2))
        self.figure.subplots_adjust(right=0.99,top=0.99)
        Canvas(pltPanel, wx.ID_ANY, self.figure) # no need to save, get this from self.figure.canvas
        self.plotaxis = self.figure.add_subplot()
        hSizer.Add(pltPanel,0,wx.EXPAND,0)
        mainSizer.Add(hSizer,1,wx.EXPAND,1)
        mainSizer.Add((-1,5))

        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        vSizer = wx.BoxSizer(wx.VERTICAL)
        self.gaugemaximum = int(maximum)
        self.gauge = wx.Gauge(self,wx.ID_ANY,range=int(self.gaugemaximum))
        self.gauge.SetValue(0)
        vSizer.Add(self.gauge,1,wx.EXPAND,1)
        vSizer.Add((-1,3))
        tSizer = wx.BoxSizer(wx.HORIZONTAL)
        tSizer.Add((-1,-1),1,wx.EXPAND,1)
        tSizer.Add(wx.StaticText(self,wx.ID_ANY,'Elapsed time:  '))
        self.elapsed = wx.StaticText(self,wx.ID_ANY,'0:00:00.0')
        tSizer.Add(self.elapsed)
        tSizer.Add((-1,-1),1,wx.EXPAND,1)
        vSizer.Add(tSizer)
        hSizer.Add(vSizer,1,wx.EXPAND,1)
        hSizer.Add((10,-1))
        btn = wx.Button(self, wx.ID_CLOSE,"Abort refinement") 
        btn.Bind(wx.EVT_BUTTON,self._onClose)
        hSizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(hSizer,0,wx.EXPAND,5)

        self.Labels.label[2].SetLabel('Start')
        self.cycleLbl = self.Labels.label[3]
        self.cycleLbl.SetLabel('Cycle ?')
        if self.trialRw:
            self.Labels.label[4].SetLabel('trial parms')
        
        self.Show()
        self.Layout()
        self.SetSizer(mainSizer)
        mainSizer.Fit(self)
        self.CenterOnParent()
        self.SendSizeEvent()
        
        self.tblCols = {}
        self.tblLbls = {}
        self.fitVals = {}
        self.trialVals = {} # used when self.trialRw is True
        self.trialcount = 0
        self.cols = 0
        self.maxCycle = 10
        self.curHist = None # number of current histogram (may be <0 for overall)
        self.seqHist = None # number of current sequential histogram (never <0)
        self.prevSeqHist = [] # previous sequential histograms that are still shown
        self.plotted = []
        self.histOff = {}
        
    def _onClose(self,event):
        '''Respond to abort button or close of window
        '''
        self.abortStatus = True
        self.Show(False)

    def Destroy(self):
        '''Destroy the window, but allow events to clear before doing so
        '''
        wx.CallAfter(wx.Dialog.Destroy,self)

    def _makeLabeledTable(self):
        '''Create two grid sizers, one with row labels and one scrolled.
        Use _xferLabeledTable to make the row heights match.
        '''        
        lblSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.Labels = wx.GridBagSizer(2,2)
        lblSizer.Add(self.Labels)
        self.RwPanel = wx.lib.scrolledpanel.ScrolledPanel(self, wx.ID_ANY, size=(180, 130),
            style = wx.SUNKEN_BORDER)
        lblSizer.Add(self.RwPanel,1,wx.EXPAND,1)

        self.Labels.label = {}
        for i in range(self.rows):
            self.Labels.label[i] = wx.StaticText(self,wx.ID_ANY,'',style=wx.ALIGN_CENTER_VERTICAL)
            self.Labels.Add(self.Labels.label[i],(i,0))
        mainsizer = wx.BoxSizer(wx.VERTICAL)
        tblsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.gridSiz = wx.GridBagSizer(2,2)
        tblsizer.Add(self.gridSiz)
        mainsizer.Add(tblsizer)
        self.RwPanel.SetSizer(mainsizer)
        self.RwPanel.SetAutoLayout(1)
        self.RwPanel.SetupScrolling()
        return lblSizer
        
    def _xferLabeledTable(self):
        '''Matches the row sizes of the labels to the row heights in the table
        '''
        for i,h in enumerate(self.gridSiz.GetRowHeights()):
            self.Labels.label[i].SetMinSize((-1,h))
        self.Labels.Layout()
        
    def SetMaxCycle(self,value):
        '''Set the maximum number of cycles or histograms (sequential fit). 
        Used to scale settings so the gauge bar completes close to 100%.
        Ignored for sequential refinements.
        '''
        if self.SeqLen: return
        self.maxCycle = value
        
    def _AddTableColumn(self,label='',col=None):
        'add a column to the Rfactor table'
        if col is None:
            self.cols += 1
            col = self.cols
        lbls = []
        for i in range(self.rows-1):
            lbls.append(wx.StaticText(self.RwPanel,wx.ID_ANY,'',style=wx.ALIGN_CENTER))
            self.gridSiz.Add(lbls[-1],(1+i,col))
        return col,lbls
        
    def SetHistogram(self,nextHist,histLbl):
        '''Set this before beginning processing of each histogram
        '''
        if nextHist is None:
            self.curHist = None
            self.msgLine1.SetLabel(histLbl)
            return
        if self.SeqLen:
            self._SetSeqHistogram(nextHist,histLbl)
            return
        col = None
        lbl = 'Hist {}'.format(nextHist)
        if nextHist == -1:  # overall fit goes in the 1st column, if shown
            col = 0
            lbl = 'Overall'
        elif nextHist == -2:  # Restraint Chi2 contribution goes in the last col
            lbl = 'Restraints'
            
        if nextHist not in self.tblCols:
            self.tblCols[nextHist],lbls = self._AddTableColumn(lbl,col)
            self.tblLbls[nextHist] = lbls
            self.tblLbls[nextHist][0].SetLabel(lbl)
            self.fitVals[nextHist] = []
        if nextHist >= 0:
            self.msgLine1.SetLabel('Fitting '+histLbl)
        self.curHist = nextHist

    def _SetSeqHistogram(self,nextHist,histLbl):
        '''Set this before beginning processing of each histogram in a sequential fit.
        Advances the completion gauge.
        '''
        if nextHist == -1:  # overall fit is not shown
            self.curHist = nextHist
            return
        elif nextHist == -2:  # Restraint Chi2 contribution goes in the 1st col
            if nextHist not in self.tblCols:
                self.tblCols[-2],lbls = self._AddTableColumn('Restraints',0)
                self.tblLbls[-2] = lbls
                self.tblLbls[-2][0].SetLabel('Restraints')
                self.fitVals[-2] = []
        elif self.seqHist != nextHist and nextHist >= 0:
            if nextHist not in self.tblCols:
                lbl = 'Hist {}'.format(nextHist)
                self.tblCols[nextHist],lbls = self._AddTableColumn(lbl,None)
                self.tblLbls[nextHist] = lbls
                self.tblLbls[nextHist][0].SetLabel(lbl)
                if len(self.prevSeqHist) < self.seqShow:
                    self.prevSeqHist.append(nextHist)
                else:
                    del self.fitVals[self.prevSeqHist[0]]                
                    del self.prevSeqHist[0]
                    self.prevSeqHist.append(nextHist)
            self.fitVals[nextHist] = []
            if -2 in self.fitVals: self.fitVals[-2] = []
            
        if nextHist >= 0:
            self.msgLine1.SetLabel('Fitting '+histLbl)
        self.curHist = nextHist
        if nextHist >= 0 and self.seqHist != nextHist:
            self.seqHist = nextHist
            self.SeqCount += 1
            self.gauge.SetValue(int(min(self.gaugemaximum,100.*self.SeqCount/self.SeqLen)))

    def _plotBar(self,h):
        'plot a vertical bar for a histogram'
        sym,lbl = self._getPlotSetting(h)
        l = len(self.fitVals[h])
        wid = 10.
        if h < 0:
            if h == -2 and -1 in self.fitVals:
                self.histOff[h] = h/wid
            else:
                self.histOff[h] = -1/wid
        else:
            self.histOff[h] = len([i for i in self.plotted if i >= 0])/wid
        self.plotaxis.bar(np.array(range(l))+self.histOff[h],self.fitVals[h],
            width=1/wid,label=lbl,color=sym)
        self.plotted.append(h)
        
    def _getPlotSetting(self,h):
        'determines how a plot is drawn'
        if h == -1:
            sym = "m" # magenta
            lbl = 'o'
        elif h == -2:
            sym = "c" # cyan
            lbl = 'r'
        else:
            symbols = ['b','r','g']
            sym = symbols[h%3]
            lbl = str(h)
        return sym,lbl

    def _SetCycleRw(self,value):
        '''Used to process an Rwp value from the :meth:`Update` method. 
        The value will be associated with the current histogram (as set
        in :meth:`SetHistogram`). If this is the 1st supplied value for
        that histogram, the value is set and displayed as as the starting 
        Rwp. If :data`:self.trialRw` is False, the values are saved to a 
        list specific to the current histogram, and are displayed and 
        plotted. When :data`:self.trialRw` is True, the Rwp values are 
        considered trial values and are only saved and plotted when 
        :meth:`AdvanceCycle` is called. 
        '''
        if self.curHist in self.fitVals: 
            cycle = len(self.fitVals[self.curHist])
        else:
            return
        if self.maxCycle and not self.SeqLen:
            self.gauge.SetValue(int(min(self.gaugemaximum,100.*cycle/self.maxCycle)))
        self.cycleLbl.SetLabel('Cycle {:}'.format(cycle))
        if cycle == 0:
            self.tblLbls[self.curHist][1].SetLabel('{:8.3g}'.format(value))
            self.RwPanel.SetupScrolling()
            self._xferLabeledTable()
        if not self.trialRw:  # show & plot current status here 
            self.fitVals[self.curHist].append(value)
            self.tblLbls[self.curHist][2].SetLabel('{:8.3g}'.format(value))
            self.plotaxis.clear()
            self.plotted = []
            if self.SeqLen:
                for h in self.fitVals:   # loop over all histograms but not all get plotted
                    if h == -1: continue
                    if h in self.prevSeqHist:
                        self._plotBar(h)
            else:
                for h in self.fitVals:   # loop over all histograms but not all get plotted
                    if h < 0 or h == self.curHist or len(self.fitVals) < 6 or (
                        self.curHist < 0 and  h==0): # plot overall & restraints, or p to 5 histograms
                        self._plotBar(h)
            if self.plotted: self.plotaxis.legend(loc=3)
            self.figure.canvas.draw()
        else: # save as trial value
            self.trialVals[self.curHist] = value
            self.tblLbls[self.curHist][3].SetLabel('{:8.3g}'.format(value))
            if self.curHist in self.plotted:
                sym,lbl = self._getPlotSetting(self.curHist)
                c = len(self.fitVals[self.curHist]) - 1 + self.histOff[self.curHist]
                self.plotaxis.plot(c,value,'o'+sym)
                self.figure.canvas.draw()
            if (self.curHist ==-1 and -2 not in self.fitVals) or self.curHist == -2:
                self.trialcount += 1
        self.Layout() # in case sizes change
        self.RwPanel.ScrollChildIntoView(self.tblLbls[self.curHist][1])

    def AdvanceCycle(self,cycle=None):
        '''Call this directly with Levenberg-Marquardt fitting after a 
        cycle completes. 
        Plots the results.
        '''
        self.plotaxis.clear()
        self.trialcount = 0
        self.plotted = []
        for h,value in self.trialVals.items():
            if h not in self.fitVals or h not in self.tblLbls: continue
            self.fitVals[h].append(value)
            self.tblLbls[h][2].SetLabel('{:8.3g}'.format(value))
            self.tblLbls[h][3].SetLabel('')
            if self.SeqLen or h < 3 or len(self.trialVals) < 6: # plot overall & restraints, or p to 5 histograms
                self._plotBar(h)
        if self.plotted: self.plotaxis.legend(loc=3)
        self.figure.canvas.draw()

    def Update(self, value=None, newmsg=""):
        '''designed to work with calls intended for wx.ProgressDialog.Update
        the value is assumed to be the current wR value for the histogram 
        selected with SetHistogram and newmsg goes into the 2nd status line. 
        '''
        if self.curHist is not None and value != 101.:
            self._SetCycleRw(value)
        if newmsg:
            self.msgLine2.SetLabel(newmsg) 
        m,s = divmod(time.time()-self.startTime,60)
        h,m = divmod(m,60)
        self.elapsed.SetLabel('{:0d}:{:02d}:{:04.1f}'.format(int(h), int(m), s))
        wx.GetApp().Yield()
        return (not self.abortStatus, True)
    
################################################################################
class downdate(wx.Dialog):
    '''Dialog to allow a user to select a version of GSAS-II to install
    svn version
    '''
    def __init__(self,parent=None):
        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, wx.ID_ANY, 'Select Version', style=style)
        pnl = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        insver = GSASIIpath.svnGetRev(local=True)
        curver = int(GSASIIpath.svnGetRev(local=False))
        label = wx.StaticText(
            pnl,  wx.ID_ANY,
            'Select a specific GSAS-II version to install'
            )
        sizer.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer1.Add(
            wx.StaticText(pnl,  wx.ID_ANY,
                          'Currently installed version: '+str(insver)),
            0, wx.ALIGN_CENTRE|wx.ALL, 5)
        sizer.Add(sizer1)
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer1.Add(
            wx.StaticText(pnl,  wx.ID_ANY,
                          'Select GSAS-II version to install: '),
            0, wx.ALIGN_CENTRE|wx.ALL, 5)
        self.spin = wx.SpinCtrl(pnl, wx.ID_ANY,size=(150,-1))
        self.spin.SetRange(1, curver)
        self.spin.SetValue(curver)
        self.Bind(wx.EVT_SPINCTRL, self._onSpin, self.spin)
        self.Bind(wx.EVT_KILL_FOCUS, self._onSpin, self.spin)
        sizer1.Add(self.spin)
        sizer.Add(sizer1)

        line = wx.StaticLine(pnl,-1, size=(-1,3), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.EXPAND|wx.ALL, 10)

        self.text = wx.StaticText(pnl,  wx.ID_ANY, "")
        sizer.Add(self.text, 0, wx.EXPAND|wx.ALL, 5)

        line = wx.StaticLine(pnl,-1, size=(-1,3), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.EXPAND|wx.ALL, 10)
        sizer.Add(
            wx.StaticText(
                pnl,  wx.ID_ANY,
                'If "Install" is pressed, your project will be saved;\n'
                'GSAS-II will exit; The specified version will be loaded\n'
                'and GSAS-II will restart. Press "Cancel" to abort.'),
            0, wx.EXPAND|wx.ALL, 10)
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(pnl, wx.ID_OK, "Install")
        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(pnl, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        pnl.SetSizer(sizer)
        sizer.Fit(self)
        self.topsizer=sizer
        self.CenterOnParent()
        self._onSpin(None)

    def _onSpin(self,event):
        'Called to load info about the selected version in the dialog'
        if event: event.Skip()
        ver = self.spin.GetValue()
        d = GSASIIpath.svnGetLog(version=ver)
        date = d.get('date','?').split('T')[0]
        s = '(Version '+str(ver)+' created '+date
        s += ' by '+d.get('author','?')+')'
        msg = d.get('msg')
        if msg: s += '\n\nComment: '+msg
        self.text.SetLabel(s)
        self.topsizer.Fit(self)

    def getVersion(self):
        'Get the version number in the dialog'
        return self.spin.GetValue()

################################################################################
class gitVersionSelector(wx.Dialog):
    '''Dialog to allow a user to select a version of GSAS-II to install
    from a git repository
    '''
    def __init__(self,parent=None):
        import git
        self.g2repo = GSASIIpath.openGitRepo(path2GSAS2)
        self.githistory = GSASIIpath.gitHistory('hash',self.g2repo)
        # patch Feb 2024: don't allow access to versions that are too old
        # since they are hard-coded to use svn
        import datetime
        tz = self.g2repo.commit(self.githistory[0]).committed_datetime.tzinfo
        cutoff = datetime.datetime(2024,2,20,tzinfo=tz)   # 20-feb-2024
        self.githistory = [h for h in self.githistory if
                        self.g2repo.commit(h).committed_datetime > cutoff]
        # end patch 
        self.initial_commit = self.g2repo.commit('HEAD')
        self.initial_commit_info = self.docCommit(self.initial_commit)
        if parent is None:
            parent = wx.GetApp().GetMainTopWindow()
        wx.Dialog.__init__(self, parent, wx.ID_ANY, 'Select GSAS-II Version',
                            style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        sizer = wx.BoxSizer(wx.VERTICAL)
        remoteupdates,localupdates = GSASIIpath.gitCountRegressions(self.g2repo)
        self.verselection = remoteupdates
        sizer.Add((-1,10))
        #label = wx.StaticText(
        #    self,  wx.ID_ANY,
        #    'Select a specific GSAS-II version to install'
        #    )
        msg = 'This allows you to revert back to a previous GSAS-II version '
        msg += 'to compare results between versions and temporarily bypass '
        msg += 'bugs. If there is something that works better in an older '
        msg += 'GSAS-II version, be sure to test the latest version and '
        msg += 'if the problem remains, report it so that it can be fixed.'
        label = wx.StaticText(self, wx.ID_ANY, msg)
        label.Wrap(400)
        sizer.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        sizer.Add((-1,20))
        
        sizer.Add(wx.StaticText(self, wx.ID_ANY,
                                    '      Currently installed version:'))
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer1.Add((50,-1))

        initpnl = wxscroll.ScrolledPanel(self, wx.ID_ANY,size=(450,90),
            style = wx.SUNKEN_BORDER)
        ssizer = wx.BoxSizer(wx.HORIZONTAL)
        txt = wx.StaticText(initpnl,  wx.ID_ANY, self.initial_commit_info)
        txt.Wrap(435)
        ssizer.Add(txt)
        initpnl.SetSizer(ssizer)
        initpnl.SetAutoLayout(1)
        initpnl.SetupScrolling()
        sizer1.Add(initpnl)
        sizer.Add(sizer1)
        
        sizer.Add((-1,20))
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer1.Add(
            wx.StaticText(self,  wx.ID_ANY,
                          'Select how many versions to regress: '),
            0, wx.ALIGN_CENTRE|wx.ALL, 5)
        self.spin = wx.SpinCtrl(self, wx.ID_ANY,size=(150,-1))
        self.spin.SetRange(-len(self.githistory)+1, 0)
        self.spin.SetValue(-self.verselection)
        self.Bind(wx.EVT_SPINCTRL, self._onSpin, self.spin)
        self.Bind(wx.EVT_KILL_FOCUS, self._onSpin, self.spin)
        sizer1.Add(self.spin)
        sizer.Add(sizer1)
        sizer.Add((-1,20))

        sizer.Add(wx.StaticText(self, wx.ID_ANY,
                                    '      Selected version to install:'))
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer1.Add((50,-1))
        
        self.spanel = wxscroll.ScrolledPanel(self, wx.ID_ANY,size=(450,90),
            style = wx.SUNKEN_BORDER)
        ssizer = wx.BoxSizer(wx.HORIZONTAL)
        self.text = wx.StaticText(self.spanel, wx.ID_ANY, "")
        ssizer.Add(self.text)
        self.spanel.SetSizer(ssizer)
        self.spanel.SetAutoLayout(1)
        self.spanel.SetupScrolling()
        sizer1.Add(self.spanel)
        sizer.Add(sizer1)

        sizer.Add((-1,20))
        sizer.Add(
            wx.StaticText(
                self,  wx.ID_ANY,
                'Press "Continue" after selecting a version to continue;\n'
                'Press "Cancel" to stop regression.'),
            0, wx.EXPAND|wx.ALL, 10)
        sizer.Add((-1,5))
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(self, wx.ID_OK, "Continue")
        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        self.SetSizer(sizer)
        sizer.Fit(self)
        siz = self.GetSize()
        siz[0] = max(siz[0],500)
        self.SetSize(siz)
        self.CenterOnParent()
        self._onSpin(None)

    def _onSpin(self,event):
        'Called to load info about the selected version in the dialog'
        if event: event.Skip()
        try:
            commit_info = self.docCommit(self.githistory[-self.spin.GetValue()])
            self.text.SetLabel(commit_info)
            self.text.Wrap(435)
        except IndexError:
            return
        self.spanel.SetAutoLayout(1)
        self.spanel.SetupScrolling()

    def getVersion(self):
        '''Gets the selected version that should be installed

        :return: returns one of three values:

         * 0: if the newest version is selected, so that the 
           installation should be updated rather than regressed
         * None: if the currently installed version is selected,
           so that nothing need be done
         * A hexsha string: the regressed version that should be 
           selected.
        '''
        if self.spin.GetValue() == 0:
            return 0
        commit = self.githistory[-self.spin.GetValue()]
        if self.g2repo.commit(commit) == self.initial_commit:
            return None
        return commit
    
    def docCommit(self,commit):
        '''Provides a string with information about a specific git commit.

        :returns: a multi-line string
        '''
        import datetime
        fmtdate = lambda c:"{:%d-%b-%Y %H:%M}".format(c.committed_datetime)
        commit = self.g2repo.commit(commit)  # converts a hash, if supplied
        msg = f'git {commit.hexsha[:10]} from {fmtdate(commit)}'
        tags = self.g2repo.git.tag('--points-at',commit).split('\n')
        if tags != ['']:
            msg += f"\ntags: {', '.join(tags)}"
        msg += '\ncomment: ' + commit.message
        return msg
    

################################################################################
class SortableLstCtrl(wx.Panel):
    '''Creates a read-only table with sortable columns. Sorting is done by 
    clicking on a column label. A triangle facing up or down is added to 
    indicate the column is sorted.

    To use, the header is labeled using 
    :meth:`PopulateHeader`, then :meth:`PopulateLine` is called for every 
    row in table and finally :meth:`SetColWidth` is called to set the column
    widths.

    :param wx.Frame parent: parent object for control
    '''
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, wx.ID_ANY)#, style=wx.WANTS_CHARS)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.list = G2LstCtrl(self, wx.ID_ANY, style=wx.LC_REPORT| wx.BORDER_SUNKEN)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        #self.SortListItems(0, True)

    def PopulateHeader(self, header, justify):
        '''Defines the column labels

        :param list header: a list of strings with header labels
        :param list justify: a list of int values where 0 causes left justification, 
          1 causes right justification, and -1 causes centering
        '''
        info = wx.ListItem()
        info.Mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info.Image = -1
        for i,(h,j) in enumerate(zip(header, justify)):
            info.Text = h
            if j > 0:
                info.Align =  wx.LIST_FORMAT_RIGHT
            elif j < 0:
                info.Align =  wx.LIST_FORMAT_CENTER
            else:
                info.Align = 0
            self.list.InsertColumn(i, info)
        listmix.ColumnSorterMixin.__init__(self.list, len(header))
        self.list.itemDataMap = {}

    def PopulateLine(self, key, data):
        '''Enters each row into the table

        :param int key: a unique int value for each line, probably should
          be sequential
        :param list data: a list of strings for each column in that row
        '''
        index = self.list.InsertItem(self.list.GetItemCount(), data[0])
        for i,d in enumerate(data[1:]):
            self.list.SetItem(index, i+1, d)
        self.list.SetItemData(index, key)
        self.list.itemDataMap[key] = data

    def SetColWidth(self,col,width=None,auto=True,minwidth=0,maxwidth=None):
        '''Sets the column width.

        :param int width: the column width in pixels
        :param bool auto: if True (default) and width is None (default) the
          width is set by the maximum width entry in the column
        :param int minwidth: used when auto is True, sets a minimum 
          column width
        :param int maxwidth: used when auto is True, sets a maximum 
          column width. Do not use with minwidth
        '''
        if width:
            self.list.SetColumnWidth(col, width)
        elif auto:
            self.list.SetColumnWidth(col, wx.LIST_AUTOSIZE)
            if minwidth and self.list.GetColumnWidth(col) < minwidth:
                self.list.SetColumnWidth(col, minwidth)
            elif maxwidth and self.list.GetColumnWidth(col) > maxwidth:
                self.list.SetColumnWidth(col, maxwidth)
        else:
            print('Error in SetColWidth: use either auto or width')

try:
    class G2LstCtrl(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin, listmix.ColumnSorterMixin):
        '''Creates a custom ListCtrl with support for images in column labels
        '''
        def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                     size=wx.DefaultSize, style=0):
            wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
            listmix.ListCtrlAutoWidthMixin.__init__(self)
            from wx.lib.embeddedimage import PyEmbeddedImage
            # from demo/images.py
            SmallUpArrow = PyEmbeddedImage(
                b"iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABHNCSVQICAgIfAhkiAAAADxJ"
                b"REFUOI1jZGRiZqAEMFGke2gY8P/f3/9kGwDTjM8QnAaga8JlCG3CAJdt2MQxDCAUaOjyjKMp"
                b"cRAYAABS2CPsss3BWQAAAABJRU5ErkJggg==")
            SmallDnArrow = PyEmbeddedImage(
                b"iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABHNCSVQICAgIfAhkiAAAAEhJ"
                b"REFUOI1jZGRiZqAEMFGke9QABgYGBgYWdIH///7+J6SJkYmZEacLkCUJacZqAD5DsInTLhDR"
                b"bcPlKrwugGnCFy6Mo3mBAQChDgRlP4RC7wAAAABJRU5ErkJggg==")
            self.il = wx.ImageList(16, 16)
            self.UpArrow = self.il.Add(SmallUpArrow.GetBitmap())
            self.DownArrow = self.il.Add(SmallDnArrow.GetBitmap())
            self.parent=parent
            self.SetImageList(self.il, wx.IMAGE_LIST_SMALL)

        def GetListCtrl(self): # needed for sorting
            return self
        def GetSortImages(self):
            #return (self.parent.DownArrow, self.parent.UpArrow)
            return (self.DownArrow, self.UpArrow)
except TypeError:
    # avoid "duplicate base class _MockObject" error in class G2LstCtrl():
    # where listmix.ListCtrlAutoWidthMixin, listmix.ColumnSorterMixin are same
    # in docs build
    class G2LstCtrl(wx.ListCtrl):
        '''Creates a custom ListCtrl with support for images in column labels
        '''
        pass
    print('docs build kludge for G2LstCtrl')
    
#### Display Help information ################################################################################
# define some globals 
htmlPanel = None
htmlFrame = None
htmlFirstUse = True
#helpLocDict = {}  # to be implemented if we ever split gsasii.html over multiple files
path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # save location of this file
def ShowHelp(helpType,frame):
    '''Called to bring up a web page for documentation.'''
    global htmlFirstUse,htmlPanel,htmlFrame
    # no defined link to use, create a default based on key
    if helpType.lower().startswith('pwdr'):
        helplink = 'gsasII-pwdr.html#'+helpType.replace(')','').replace('(','_').replace(' ','_')
    elif helpType.lower().startswith('phase'):
        helplink = 'gsasII-phase.html#'+helpType.replace(')','').replace('(','_').replace(' ','_')
    elif helpType.lower().startswith('hist/phase'):
        helplink = 'gsasII-phase.html#Phase-Data'
    elif helpType:
        helplink = 'gsasII.html#'+helpType.replace(')','').replace('(','_').replace(' ','_')
    else: 
        helplink = 'gsasII.html'
    # determine if a web browser or the internal viewer should be used for help info
    if GSASIIpath.GetConfigValue('Help_mode'):
        helpMode = GSASIIpath.GetConfigValue('Help_mode')
    else:
        helpMode = 'browser'
    if helpMode == 'internal':
        helplink = os.path.join(path2GSAS2,'help',helplink)
        try:
            htmlPanel.LoadFile(helplink)
            htmlFrame.Raise()
        except:
            htmlFrame = wx.Frame(frame, -1, size=(610, 510))
            htmlFrame.Show(True)
            htmlFrame.SetTitle("HTML Window") # N.B. reset later in LoadFile
            htmlPanel = MyHtmlPanel(htmlFrame,-1)
            htmlPanel.LoadFile(helplink)
    else:
        import webbrowser     # postpone this until now for quicker startup
        wb = webbrowser
        if sys.platform == "darwin": # on Mac, use a OSXscript so that file anchors work
            # Get the default browser, this will fail in py2.7 and might fail, so 
            # use safari as a backup
            appleScript = '''    
    use framework "AppKit"
    use AppleScript version "2.4"
    use scripting additions

    property NSWorkspace : a reference to current application's NSWorkspace
    property NSURL : a reference to current application's NSURL

    set wurl to NSURL's URLWithString:"https://www.apple.com"
    set thisBrowser to (NSWorkspace's sharedWorkspace)'s ¬
                        URLForApplicationToOpenURL:wurl
    set appname to (thisBrowser's absoluteString)'s lastPathComponent()'s ¬
                    stringByDeletingPathExtension() as text
    return appname as text
            '''
            import subprocess
            try:
                browser = subprocess.check_output(["osascript","-e",appleScript], encoding='UTF-8').strip()
                wb = webbrowser.MacOSXOSAScript(browser)
            except:
                wb = webbrowser.MacOSXOSAScript('safari')
        # open the link
        helplink = os.path.join(path2GSAS2,'help',helplink)
        pfx = "file://"
        if sys.platform.lower().startswith('win'):
            pfx = ''
        #if GSASIIpath.GetConfigValue('debug'): print 'DBG_Help link=',pfx+helplink
        if htmlFirstUse:
            wb.open_new(pfx+helplink)
            htmlFirstUse = False
        else:
            wb.open(pfx+helplink, new=0, autoraise=True)

def ShowWebPage(URL,frame,browser=False):
    '''Called to show a tutorial web page.

    :param str URL: web page URL
    :param wx.Frame frame: parent window (or None)
    :param bool browser: If True, forces the page to be opened in a web 
      browser, regardless of the ``Help_mode`` config setting. 
    '''
    global htmlFirstUse,htmlPanel,htmlFrame
    # determine if a web browser or the internal viewer should be used for help info
    if GSASIIpath.GetConfigValue('Help_mode') and not browser:
        helpMode = GSASIIpath.GetConfigValue('Help_mode')
    else:
        helpMode = 'browser'
    if helpMode == 'internal':
        try:
            htmlPanel.LoadFile(URL)
            htmlFrame.Raise()
        except:
            htmlFrame = wx.Frame(frame, -1, size=(610, 510))
            htmlFrame.Show(True)
            htmlFrame.SetTitle("HTML Window") # N.B. reset later in LoadFile
            htmlPanel = MyHtmlPanel(htmlFrame,-1)
            htmlPanel.LoadFile(URL)
    else:
        import webbrowser     # postpone this until now for quicker startup
        if URL.startswith('http'): 
            pfx = ''
        elif sys.platform.lower().startswith('win'):
            pfx = ''
        else:
            pfx = "file://"
        if htmlFirstUse:
            webbrowser.open_new(pfx+URL)
            htmlFirstUse = False
        else:
            webbrowser.open(pfx+URL, new=0, autoraise=True)

#### Tutorials support ################################################################################
G2TutURL = "https://advancedphotonsource.github.io/GSAS-II-tutorials/"

tutorialCatalog = [l for l in tutorialIndex if len(l) >= 3]
# A catalog of GSAS-II tutorials generated from the table in :data:`tutorialIndex`
def OpenTutorial(parent):
    if GSASIIpath.HowIsG2Installed().startswith('git'):
        return OpenGitTutorial(parent)
    else:
        return OpenSvnTutorial(parent)

class OpenSvnTutorial(wx.Dialog):
    '''Open a tutorial web page, optionally copying the web page, screen images and
    data file(s) to the local disk.
    '''
    
    def __init__(self,parent):
        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, wx.ID_ANY, 'Open Tutorial', style=style)
        self.G2frame = self.frame = parent
        pnl = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)        
        label = wx.StaticText(
            pnl,  wx.ID_ANY,
            'Select the tutorial to be run and the mode of access'
            )
        msg = '''GSAS-II tutorials and their sample data files
        require a fair amount of storage space; few users will
        use all of them. This dialog allows users to load selected
        tutorials (along with their sample data) to their computer;
        optionally all tutorials can be downloaded.

        Downloaded tutorials can be viewed and run without internet
        access. Tutorials can also be viewed without download, but
        users will need to download the sample data files manually.

        The location used to download tutorials is set using the
        "Set download location" which is saved as the "Tutorial_location"
        configuration option see File/Preference or the
        config_example.py file. System managers can select to have
        tutorial files installed at a shared location.
        '''
        self.SetTutorialPath()
        hlp = HelpButton(pnl,msg)
        sizer1.Add((-1,-1),1, wx.EXPAND, 0)
        sizer1.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 0)
        sizer1.Add((-1,-1),1, wx.EXPAND, 0)
        sizer1.Add(hlp)
        sizer.Add(sizer1,0,wx.EXPAND|wx.ALL,0)
        sizer.Add((10,10))
        sizer0 = wx.BoxSizer(wx.HORIZONTAL)        
        sizer1a = wx.BoxSizer(wx.VERTICAL)
        sizer1b = wx.BoxSizer(wx.VERTICAL)
        btn = wx.Button(pnl, wx.ID_ANY, "Download a tutorial and view")
        btn.Bind(wx.EVT_BUTTON, self.SelectAndDownload)
        sizer1a.Add(btn,0)
        btn = wx.Button(pnl, wx.ID_ANY, "Select from downloaded tutorials")
        btn.Bind(wx.EVT_BUTTON, self.onSelectDownloaded)
        sizer1a.Add(btn,0)
        btn = wx.Button(pnl, wx.ID_ANY, "Browse tutorial on web")
        btn.Bind(wx.EVT_BUTTON, self.onWebBrowse)
        sizer1a.Add(btn,0)
        btn = wx.Button(pnl, wx.ID_ANY, "Update downloaded tutorials")
        btn.Bind(wx.EVT_BUTTON, self.UpdateDownloaded)
        sizer1b.Add(btn,0)
        btn = wx.Button(pnl, wx.ID_ANY, "Download all tutorials")
        btn.Bind(wx.EVT_BUTTON, self.DownloadAll)
        sizer1b.Add(btn,0)
        sizer0.Add(sizer1a,0,wx.EXPAND|wx.ALL,0)
        sizer0.Add(sizer1b,0,wx.EXPAND|wx.ALL,0)
        sizer.Add(sizer0,5,wx.EXPAND|wx.ALL,5)
        
        sizer.Add((10,10))
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(pnl, wx.ID_ANY, "Set download location")
        btn.Bind(wx.EVT_BUTTON, self.SelectDownloadLoc)
        sizer1.Add(btn,0,WACV)
        self.dataLoc = wx.StaticText(pnl, wx.ID_ANY,self.tutorialPath)
        sizer1.Add(self.dataLoc,0,WACV)
        sizer.Add(sizer1)
        
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(pnl, wx.ID_CANCEL,"Done")
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        pnl.SetSizer(sizer)
        sizer.Fit(self)
        self.topsizer=sizer
        self.CenterOnParent()

    def SetTutorialPath(self):
        '''Get the tutorial location if set; if not pick a default
        directory in a logical place
        '''
        # has the user set a location and is it valid?
        if GSASIIpath.GetConfigValue('Tutorial_location'):
            tutorialPath = os.path.abspath(GSASIIpath.GetConfigValue('Tutorial_location'))
            if os.path.exists(tutorialPath):
                self.tutorialPath = tutorialPath
                return
            try:
                os.makedirs(tutorialPath)
                if os.path.exists(tutorialPath):
                    self.tutorialPath = tutorialPath
                    return
            except:
                print('Unable to use Tutorial_location config setting',
                          tutorialPath)
        # try a system-specific location
        if (sys.platform.lower().startswith('win')):
            for p in ('Documents','My Documents',''):
                if os.path.exists(os.path.abspath(os.path.expanduser(
                      os.path.join('~',p)))):
                    tutorialPath = os.path.abspath(os.path.expanduser(
                      os.path.join('~',p,'G2tutorials')))
                    if os.path.exists(tutorialPath):
                        self.tutorialPath = tutorialPath
                        return
                    try:
                        os.makedirs(tutorialPath)
                        if os.path.exists(tutorialPath):
                            self.tutorialPath = tutorialPath
                            return
                    except:
                        pass
        else:
            tutorialPath = os.path.abspath(os.path.expanduser(
                    os.path.join('~','G2tutorials')))
            if os.path.exists(tutorialPath):
                self.tutorialPath = tutorialPath
                return
            try:
                os.makedirs(tutorialPath)
                if os.path.exists(tutorialPath):
                    self.tutorialPath = tutorialPath
                    return
            except:
                pass
        # no success so far, use current working directory
        tutorialPath = os.path.abspath(os.path.join(os.getcwd(),'G2tutorials'))
        if os.path.exists(tutorialPath):
            self.tutorialPath = tutorialPath
            return
        try:
            os.makedirs(tutorialPath)
            if os.path.exists(tutorialPath):
                self.tutorialPath = tutorialPath
                return
        except:
            pass
        # nothing worked, set self.tutorialPath with os.getcwd() and hope for the best
        print('Warning: Unable to set a TutorialPath, using',os.getcwd())
        tutorialPath = os.getcwd()

    def SelectAndDownload(self,event):
        '''Make a list of all tutorials on web and allow user to choose one to
        download and then view
        '''
        indices = [j for j,i in enumerate(tutorialCatalog)
            if not os.path.exists(os.path.join(self.tutorialPath,i[0],i[1]))]
        if not indices:
            G2MessageBox(self,'All tutorials are downloaded','None to download')
            return
#        choices = [tutorialCatalog[i][2] for i in indices]
#        selected = self.ChooseTutorial(choices)
        choices2 = [tutorialCatalog[i][2:4] for i in indices]
        selected = self.ChooseTutorial2(choices2)
        if selected is None: return
        j = indices[selected]
        fullpath = os.path.join(self.tutorialPath,tutorialCatalog[j][0],tutorialCatalog[j][1])
        fulldir = os.path.join(self.tutorialPath,tutorialCatalog[j][0])
        #URL = G2BaseURL+'/Tutorials/'+tutorialCatalog[j][0]+'/'
        URL = G2TutURL + tutorialCatalog[j][0]+'/'
        if GSASIIpath.svnInstallDir(URL,fulldir):
            ShowWebPage(fullpath,self.frame)
        else:
            G2MessageBox(self,'Error downloading tutorial','Download error')
        self.EndModal(wx.ID_OK)
        self.G2frame.TutorialImportDir = os.path.join(self.tutorialPath,tutorialCatalog[j][0],'data')

    def onSelectDownloaded(self,event):
        '''Select a previously downloaded tutorial
        '''
        indices = [j for j,i in enumerate(tutorialCatalog)
            if os.path.exists(os.path.join(self.tutorialPath,i[0],i[1]))]
        if not indices:
            G2MessageBox(self,
                         'There are no downloaded tutorials in '+self.tutorialPath,
                         'None downloaded')
            return
#        choices = [tutorialCatalog[i][2] for i in indices]
#        selected = self.ChooseTutorial(choices)
        choices2 = [tutorialCatalog[i][2:4] for i in indices]
        selected = self.ChooseTutorial2(choices2)
        if selected is None: return
        j = indices[selected]
        fullpath = os.path.join(self.tutorialPath,tutorialCatalog[j][0],tutorialCatalog[j][1])
        self.EndModal(wx.ID_OK)
        ShowWebPage(fullpath,self.frame)
        self.G2frame.TutorialImportDir = os.path.join(self.tutorialPath,tutorialCatalog[j][0],'data')
        
    def onWebBrowse(self,event):
        '''Make a list of all tutorials on web and allow user to view one.
        '''
#        choices = [i[2] for i in tutorialCatalog]
#        selected = self.ChooseTutorial(choices)
        choices2 = [i[2:4] for i in tutorialCatalog]
        selected = self.ChooseTutorial2(choices2)
        if selected is None: return        
        tutdir = tutorialCatalog[selected][0]
        tutfil = tutorialCatalog[selected][1]
        # open web page remotely, don't worry about data
        #URL = G2BaseURL+'/Tutorials/'+tutdir+'/'+tutfil
        URL = G2TutURL + tutdir + '/' +tutfil
        self.EndModal(wx.ID_OK)
        ShowWebPage(URL,self.frame)
        
    def ChooseTutorial2(self,choices):
        '''Select tutorials from a two-column table, when possible
        '''
        lbls = ('tutorial name (indent indicates previous is required)','description')
        colWidths=[400,400]
        dlg = MultiColumnSelection(self,'select tutorial',lbls,choices,colWidths)
        selection = dlg.Selection
        dlg.Destroy()
        if selection is not None: # wx from EPD Python
            if selection == -1: return
            return selection
        else:
            return self.ChooseTutorial([i[0] for i in choices])
        
    def ChooseTutorial(self,choices):
        '''choose a tutorial from a list
        (will eventually only be used with very old wxPython
        '''
        def onDoubleClick(event):
            'double-click closes the dialog'
            dlg.EndModal(wx.ID_OK)
        dlg = wx.Dialog(self,wx.ID_ANY,
                        'Select a tutorial to view. NB: indented ones require prerequisite',
                        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        pnl = wx.Panel(dlg)
        sizer = wx.BoxSizer(wx.VERTICAL)
        listbox = wx.ListBox(pnl, wx.ID_ANY, choices=choices,size=(450, 100),style=wx.LB_SINGLE)
        sizer.Add(listbox,1,wx.EXPAND|wx.ALL,1)
        listbox.Bind(wx.EVT_LISTBOX_DCLICK, onDoubleClick)
        sizer.Add((10,10))
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(pnl, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        OKbtn = wx.Button(pnl, wx.ID_OK)
        OKbtn.SetDefault()
        btnsizer.AddButton(OKbtn)
        btnsizer.Realize()
        sizer.Add((-1,5))
        sizer.Add(btnsizer,0,wx.ALIGN_RIGHT,50)
        
        pnl.SetSizer(sizer)
        sizer.Fit(dlg)
        dlg.CenterOnParent()
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        selected = listbox.GetSelection()
        dlg.Destroy()
        wx.Yield() # close window right away so user sees something happen
        if selected < 0: return
        return selected

    def UpdateDownloaded(self,event):
        'Find the downloaded tutorials and run an svn update on them'
        updated = 0
        wx.BeginBusyCursor()
        for i in tutorialCatalog:
            if not os.path.exists(os.path.join(self.tutorialPath,i[0],i[1])): continue
            print('Updating '+i[0])
            GSASIIpath.svnUpdateDir(os.path.join(self.tutorialPath,i[0]))
            updated += 1
        wx.EndBusyCursor()
        if not updated:
            G2MessageBox(self,'Warning, you have no downloaded tutorials','None Downloaded')
        else:
            G2MessageBox(self,'{} updates completed'.format(updated),'Updates done')
        #self.EndModal(wx.ID_OK)
        
    def DownloadAll(self,event):
        'Download or update all tutorials'
        fail = ''
        for i in tutorialCatalog:
            if os.path.exists(os.path.join(self.tutorialPath,i[0],i[1])):
                print('Updating '+i[0])
                GSASIIpath.svnUpdateDir(os.path.join(self.tutorialPath,i[0]))
            else:
                fulldir = os.path.join(self.tutorialPath,i[0])
                #URL = G2BaseURL+'/Tutorials/'+i[0]+'/'
                URL = G2TutURL + i[0] + '/'
                if not GSASIIpath.svnInstallDir(URL,fulldir):
                    if fail: fail += ', '
                    fail += i[0]
        if fail: 
            G2MessageBox(self,'Error downloading tutorial(s)\n\t'+fail,'Download error')
        self.EndModal(wx.ID_OK)
                    
    def SelectDownloadLoc(self,event):
        '''Select a download location,
        Cancel resets to the default
        '''
        dlg = wx.DirDialog(self, "Choose a directory for tutorial downloads:",
            defaultPath=self.tutorialPath)#,style=wx.DD_DEFAULT_STYLE)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            pth = dlg.GetPath()
        finally:
            dlg.Destroy()

        if not os.path.exists(pth):
            try:
                os.makedirs(pth)    #failing for no obvious reason
            except OSError:
                msg = 'The selected directory is not valid.\n\t'
                msg += pth
                msg += '\n\nAn attempt to create the directory failed'
                G2MessageBox(self.frame,msg)
                return
        if os.path.exists(os.path.join(pth,"help")) and os.path.exists(os.path.join(pth,"Exercises")):
            print("Note that you may have old tutorial files in the following directories")
            print('\t'+os.path.join(pth,"help"))
            print('\t'+os.path.join(pth,"Exercises"))
            print('Subdirectories in the above can be deleted to save space\n\n')
        self.tutorialPath = pth
        self.dataLoc.SetLabel(self.tutorialPath)
        if GSASIIpath.GetConfigValue('Tutorial_location') == pth: return
        vars = GetConfigValsDocs()
        try:
            vars['Tutorial_location'][1] = pth
            if GSASIIpath.GetConfigValue('debug'): print('DBG_Saving Tutorial_location: '+pth)
            GSASIIpath.SetConfigValue(vars)
            SaveConfigVars(vars)
        except KeyError:
            pass
        
class OpenGitTutorial(wx.Dialog):
    '''Open a tutorial web page from the git repository, 
    optionally copying the tutorial's exercise data file(s) to 
    the local disk.
    '''
    
    def __init__(self,parent):
        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, wx.ID_ANY, 'Open Tutorial', style=style)
        self.G2frame = self.frame = parent
        pnl = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        msg = ('To open a tutorial, select below the mode to be used. '+
               'With either choice you can select a tutorial and view it '+
               'in a web browser, but the 2nd option offers the option '+
               'to automatically download the data files needed to run '+
               'that tutorial.')
        label = wx.StaticText(pnl,  wx.ID_ANY, msg)
        label.Wrap(450)
        msg = '''The data files needed to run the GSAS-II tutorials
        require a fair amount of storage space; few users will
        use all of them. This dialog allows you to open a
        tutorial in a web browser and select if you want the data 
        needed to run that exercise to be downloaded to your computer.

        The location used to download tutorials is set using the
        "Set download location" which is saved as the "Tutorial_location"
        configuration option see File/Preference or the
        config_example.py file.
        '''
        self.SetTutorialPath()
        hlp = HelpButton(pnl,msg)
        sizer1.Add((-1,-1),1, wx.EXPAND, 0)
        sizer1.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        sizer1.Add((-1,-1),1, wx.EXPAND, 0)
        sizer1.Add(hlp)
        sizer.Add(sizer1,0,wx.EXPAND|wx.ALL,0)
        sizer.Add((10,10))
        sizer0 = wx.BoxSizer(wx.HORIZONTAL)        
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        btn = wx.Button(pnl, wx.ID_ANY, "Select a tutorial to view")
        btn.Bind(wx.EVT_BUTTON, self.onWebBrowse)
        sizer1.Add(btn,0)
        btn = wx.Button(pnl, wx.ID_ANY, "Select a tutorial to view and download its data files")
        btn.Bind(wx.EVT_BUTTON, self.SelectAndDownload)
        sizer1.Add(btn,0,wx.TOP,5)
        sizer0.Add(sizer1,0,wx.EXPAND|wx.ALL,0)
        sizer.Add(sizer0,5,wx.EXPAND|wx.ALL,5)
        
        sizer.Add((10,10))
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(pnl, wx.ID_ANY, "Set download location")
        btn.Bind(wx.EVT_BUTTON, self.SelectDownloadLoc)
        sizer1.Add(btn,0,WACV|wx.RIGHT|wx.LEFT,5)
        self.dataLoc = wx.StaticText(pnl, wx.ID_ANY,self.tutorialPath)
        sizer1.Add(self.dataLoc,0,WACV)
        sizer.Add(sizer1)
        
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(pnl, wx.ID_CANCEL,"Close")
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        pnl.SetSizer(sizer)
        sizer.Fit(self)
        self.topsizer=sizer
        self.CenterOnParent()

    def SetTutorialPath(self):
        '''Get the tutorial location if set; if not pick a default
        directory in a logical place
        '''
        # has the user set a location and is it valid?
        if GSASIIpath.GetConfigValue('Tutorial_location'):
            tutorialPath = os.path.abspath(GSASIIpath.GetConfigValue('Tutorial_location'))
            if os.path.exists(tutorialPath):
                self.tutorialPath = tutorialPath
                return
            try:
                os.makedirs(tutorialPath)
                if os.path.exists(tutorialPath):
                    self.tutorialPath = tutorialPath
                    return
            except:
                print('Unable to use Tutorial_location config setting',
                          tutorialPath)
        # try a system-specific location
        if (sys.platform.lower().startswith('win')):
            for p in ('Documents','My Documents',''):
                if os.path.exists(os.path.abspath(os.path.expanduser(
                      os.path.join('~',p)))):
                    tutorialPath = os.path.abspath(os.path.expanduser(
                      os.path.join('~',p,'G2tutorials')))
                    if os.path.exists(tutorialPath):
                        self.tutorialPath = tutorialPath
                        return
                    try:
                        os.makedirs(tutorialPath)
                        if os.path.exists(tutorialPath):
                            self.tutorialPath = tutorialPath
                            return
                    except:
                        pass
        else:
            tutorialPath = os.path.abspath(os.path.expanduser(
                    os.path.join('~','G2tutorials')))
            if os.path.exists(tutorialPath):
                self.tutorialPath = tutorialPath
                return
            try:
                os.makedirs(tutorialPath)
                if os.path.exists(tutorialPath):
                    self.tutorialPath = tutorialPath
                    return
            except:
                pass
        # no success so far, use current working directory
        tutorialPath = os.path.abspath(os.path.join(os.getcwd(),'G2tutorials'))
        if os.path.exists(tutorialPath):
            self.tutorialPath = tutorialPath
            return
        try:
            os.makedirs(tutorialPath)
            if os.path.exists(tutorialPath):
                self.tutorialPath = tutorialPath
                return
        except:
            pass
        # nothing worked, set self.tutorialPath with os.getcwd() and hope for the best
        print('Warning: Unable to set a TutorialPath, using',os.getcwd())
        tutorialPath = os.getcwd()

    def SelectAndDownload(self,event):
        '''Shows a list of all tutorials so user can select one to view.
        The data files associated with that directory are then downloaded.
        '''
        tutdir = self.onWebBrowse(event)
        GSASIIpath.downloadDirContents([tutdir,'data'],self.tutorialPath)

    def onWebBrowse(self,event):
        '''Shows a list of all tutorials so user can select one to view.

        :returns: the name of the directory where the tutorial is located,
          which is used if called from :meth:`SelectAndDownload`.
        '''
        choices2 = [i[2:4] for i in tutorialCatalog]
        selected = self.ChooseTutorial2(choices2)
        if selected is None: return        
        tutdir = tutorialCatalog[selected][0]
        tutfil = tutorialCatalog[selected][1]
        # open web page remotely, don't worry about data
        #URL = G2BaseURL+'/Tutorials/'+tutdir+'/'+tutfil
        URL = G2TutURL + tutdir + '/' + tutfil
        wx.CallAfter(self.EndModal,wx.ID_OK)
        ShowWebPage(URL,self.frame)
        return tutdir
        
    def ChooseTutorial2(self,choices):
        '''Select tutorials from a two-column table, when possible
        '''
        lbls = ('tutorial name (indent indicates previous is required)','description')
        colWidths=[400,400]
        dlg = MultiColumnSelection(self,'select tutorial',lbls,choices,colWidths)
        selection = dlg.Selection
        dlg.Destroy()
        if selection is not None:
            if selection == -1: return
            return selection

    def SelectDownloadLoc(self,event):
        '''Select a download location,
        Cancel resets to the default
        '''
        dlg = wx.DirDialog(self, "Choose a directory for tutorial downloads:",
            defaultPath=self.tutorialPath)#,style=wx.DD_DEFAULT_STYLE)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            pth = dlg.GetPath()
        finally:
            dlg.Destroy()

        if not os.path.exists(pth):
            try:
                os.makedirs(pth)    #failing for no obvious reason
            except OSError:
                msg = 'The selected directory is not valid.\n\t'
                msg += pth
                msg += '\n\nAn attempt to create the directory failed'
                G2MessageBox(self.frame,msg)
                return
        if os.path.exists(os.path.join(pth,"help")) and os.path.exists(os.path.join(pth,"Exercises")):
            print("Note that you may have old tutorial files in the following directories")
            print('\t'+os.path.join(pth,"help"))
            print('\t'+os.path.join(pth,"Exercises"))
            print('Subdirectories in the above can be deleted to save space\n\n')
        self.tutorialPath = pth
        self.dataLoc.SetLabel(self.tutorialPath)
        if GSASIIpath.GetConfigValue('Tutorial_location') == pth: return
        vars = GetConfigValsDocs()
        try:
            vars['Tutorial_location'][1] = pth
            if GSASIIpath.GetConfigValue('debug'): print('DBG_Saving Tutorial_location: '+pth)
            GSASIIpath.SetConfigValue(vars)
            SaveConfigVars(vars)
        except KeyError:
            pass

### Autoload PWDR files ################################################################################
AutoLoadWindow = None

def AutoLoadFiles(G2frame,FileTyp='pwd'):
    def OnBrowse(event):
        '''Responds when the Browse button is pressed to load a file.
        The routine determines which button was pressed and gets the
        appropriate file type and loads it into the appropriate place
        in the dict.
        '''
        if btn3 == event.GetEventObject():
            d = wx.DirDialog(dlg,
                    'Select directory for input files',
                    Settings['indir'],wx.DD_DEFAULT_STYLE)
            d.CenterOnParent()
            try:
                if d.ShowModal() == wx.ID_OK:
                    Settings['indir'] = d.GetPath()
                    fInp3.SetValue(Settings['indir'])
            finally:
                d.Destroy()
        elif btn4 == event.GetEventObject():
            extList = 'GSAS iparm file (*.prm,*.inst,*.ins,.instprm)|*.prm;*.inst;*.ins;*.instprm'
            d = wx.FileDialog(dlg,
                'Choose instrument parameter file',
                '', '',extList, wx.FD_OPEN)
            if os.path.exists(Settings['instfile']):
                dr,f = os.path.split(Settings['instfile'])
                d.SetDirectory(dr)
                d.SetFilename(f)
            try:
                if d.ShowModal() == wx.ID_OK:
                    Settings['instfile'] = d.GetPath()
                    fInp4.SetValue(Settings['instfile'])
                    # change the "read from" directory if defaulted
                    if Settings['indir'] == os.getcwd():
                        Settings['indir'] = os.path.dirname(d.GetPath())
                        fInp3.SetValue(Settings['indir'])
            finally:
                d.Destroy()
        TestInput()
        
    def OnFileOfFiles(event):
        '''Read from a list of files and add those files in the order
        specified in that file.
        '''
        # get a list of existing histograms
        if FileTyp == 'pwd':
            treePrfx = 'PWDR '
        else:
            treePrfx = 'PDF  '
        ReadList = []
        if G2frame.GPXtree.GetCount():
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                name = G2frame.GPXtree.GetItemText(item)
                if name.startswith(treePrfx) and name not in ReadList:
                    ReadList.append(name)
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
        Settings['ReadList'] = ReadList
        extList = 'text file (*.txt,*.csv)|*.txt;*.csv'
        d = wx.FileDialog(dlg,
                'Choose a text file with file names',
                '', '',extList, wx.FD_OPEN)
        filelist = []
        try:
            if d.ShowModal() == wx.ID_OK:
                if not os.path.exists(Settings['instfile']): return
                if not os.path.exists(d.GetPath()): return
                with open(d.GetPath(),'r') as fp:
                    for line in fp: # split lines at comma/tab, strip quotes, etc
                        if line.startswith('#'): continue
                        line = line.split(',')[0]
                        line = line.split('\t')[0]
                        f = line.replace('"','').replace("'",'').strip()
                        if not os.path.exists(f) and not os.path.abspath(f):
                            f = os.path.join(Settings['indir'],f)
                        if not os.path.exists(f):
                            print(f'Skipping file {f}, not found')
                        else:
                            filelist.append(f)
                G2frame.CheckNotebook()
                RunTimerPWDR(None,filelist)
                wx.CallAfter(dlg.Destroy)
        finally:
            d.Destroy()
    def onSetFmtSelection():
        extSel.Clear()
        extSel.AppendItems(fileReaders[Settings['fmt']].extensionlist)
        Settings['extStr'] = fileReaders[Settings['fmt']].extensionlist[0]
        extSel.SetSelection(0)
        onSetExtSelection()
    def onSetExtSelection():
        Settings['filter'] = os.path.splitext(Settings['filter'])[0] + Settings['extStr']
        flterInp.SetValue(Settings['filter'])
        TestInput()
    def OnQuit(event):
        Settings['timer'].Stop()
        wx.CallAfter(dlg.Destroy)
    def TestInput(*args,**kwargs):
        valid = True
        if not os.path.exists(Settings['indir']):
            valid = False
        if FileTyp == 'pwd' and not os.path.exists(Settings['instfile']):
            valid = False
        btnstart.Enable(valid)
        FofFbtn.Enable(valid)
    def OnStart(event):
        if btnstart.GetLabel() == 'Pause':
            Settings['timer'].Stop()
            btnstart.SetLabel('Continue')
            return
        else:
            btnstart.SetLabel('Pause')
        if Settings['timer'].IsRunning(): return
        PollTime = 1 # sec
        G2frame.CheckNotebook()
        Settings['timer'].Start(int(1000*PollTime),oneShot=False)

        # get a list of existing histograms
        if FileTyp == 'pwd':
            treePrfx = 'PWDR '
        else:
            treePrfx = 'PDF  '
        ReadList = []
        if G2frame.GPXtree.GetCount():
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                name = G2frame.GPXtree.GetItemText(item)
                if name.startswith(treePrfx) and name not in ReadList:
                    ReadList.append(name)
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
        Settings['ReadList'] = ReadList
    def RunTimerPWDR(event,filelist=None):
        if filelist is None:
            if GSASIIpath.GetConfigValue('debug'):
                import datetime
                print ("DBG_Timer tick at {:%d %b %Y %H:%M:%S}\n".format(datetime.datetime.now()))
            filelist = glob.glob(os.path.join(Settings['indir'],Settings['filter']))
            if not filelist: return
        #if GSASIIpath.GetConfigValue('debug'): print(filelist)
        Id = None
        for f in filelist:
            if f in Settings['filesread']: continue
            Settings['filesread'].append(f)
            rd = fileReaders[Settings['fmt']]
            rd.ReInitialize()
            if not rd.ContentsValidator(f):
                Settings['timer'].Stop()
                btnstart.SetLabel('Continue')
                G2MessageBox(dlg,'Error in reading file {}: {}'.format(
                    f, rd.errors))
                return
            #if len(rd.selections) > 1:
            #    G2fil.G2Print('Warning: Skipping file {}: multibank not yet implemented'.format(f))
            #    continue
            block = 0
            rdbuffer = {}
            repeat = True
            while repeat:
                repeat = False
                try:
                    flag = rd.Reader(f,buffer=rdbuffer, blocknum=block)
                except:
                    flag = False
                if flag:
                    rd.readfilename = f
                    if rd.warnings:
                        G2fil.G2Print("Read warning by", rd.formatName,
                                          "reader:",
                                          rd.warnings)
                    elif not block:
                        G2fil.G2Print("{} read by Reader {}"
                              .format(f,rd.formatName))
                    else:
                        G2fil.G2Print("{} block # {} read by Reader {}"
                                .format(f,block,rd.formatName))
                    block += 1    
                    repeat = rd.repeat
                else:
                    G2fil.G2Print("Warning: {} Reader failed to read {}"
                                      .format(rd.formatName,f))
                Iparm1, Iparm2 = G2sc.load_iprms(Settings['instfile'],rd)    
                if 'phoenix' in wx.version():
                    HistName = 'PWDR '+rd.idstring
                else:
                    HistName = 'PWDR '+G2obj.StripUnicode(rd.idstring,'_')
                # make new histogram names unique
                HistName = G2obj.MakeUniqueLabel(HistName,Settings['ReadList'])
                Settings['ReadList'].append(HistName)
                # put into tree
                Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=HistName)
                if 'T' in Iparm1['Type'][0]:
                    if not rd.clockWd and rd.GSAS:
                        rd.powderdata[0] *= 100.        #put back the CW centideg correction
                    cw = np.diff(rd.powderdata[0])
                    rd.powderdata[0] = rd.powderdata[0][:-1]+cw/2.
                    if rd.GSAS:     #NB: old GSAS wanted intensities*CW even if normalized!
                        npts = min(len(rd.powderdata[0]),len(rd.powderdata[1]),len(cw))
                        rd.powderdata[1] = rd.powderdata[1][:npts]/cw[:npts]
                        rd.powderdata[2] = rd.powderdata[2][:npts]*cw[:npts]**2  #1/var=w at this point
                    else:       #NB: from topas/fullprof type files
                        rd.powderdata[1] = rd.powderdata[1][:-1]
                        rd.powderdata[2] = rd.powderdata[2][:-1]
                    if 'Itype' in Iparm2:
                        Ibeg = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][0])
                        Ifin = np.searchsorted(rd.powderdata[0],Iparm2['Tminmax'][1])
                        rd.powderdata[0] = rd.powderdata[0][Ibeg:Ifin]
                        YI,WYI = G2pwd.calcIncident(Iparm2,rd.powderdata[0])
                        rd.powderdata[1] = rd.powderdata[1][Ibeg:Ifin]/YI
                        var = 1./rd.powderdata[2][Ibeg:Ifin]
                        var += WYI*rd.powderdata[1]**2
                        var /= YI**2
                        rd.powderdata[2] = 1./var
                    rd.powderdata[3] = np.zeros_like(rd.powderdata[0])
                    rd.powderdata[4] = np.zeros_like(rd.powderdata[0])
                    rd.powderdata[5] = np.zeros_like(rd.powderdata[0])
                Ymin = np.min(rd.powderdata[1])                 
                Ymax = np.max(rd.powderdata[1])                 
                valuesdict = {
                    'wtFactor':1.0,
                    'Dummy':False,
                    'ranId':ran.randint(0,sys.maxsize),
                    'Offset':[0.0,0.0],'delOffset':0.02*Ymax,'refOffset':-.1*Ymax,'refDelt':0.1*Ymax,
                    'Yminmax':[Ymin,Ymax]
                    }
                # apply user-supplied corrections to powder data
                if 'CorrectionCode' in Iparm1:
                    print('Warning: CorrectionCode from instprm file not applied')
                rd.Sample['ranId'] = valuesdict['ranId'] # this should be removed someday
                G2frame.GPXtree.SetItemPyData(Id,[valuesdict,rd.powderdata])
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Comments'),
                    rd.comments)
                Tmin = min(rd.powderdata[0])
                Tmax = max(rd.powderdata[0])
                Tmin1 = Tmin
                if 'NT' in Iparm1['Type'][0] and G2lat.Pos2dsp(Iparm1,Tmin) < 0.4:                
                    Tmin1 = G2lat.Dsp2pos(Iparm1,0.4)
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Limits'),
                    rd.pwdparms.get('Limits',[(Tmin,Tmax),[Tmin1,Tmax]])
                    )
                G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame,Id,'Limits')
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Background'),
                    rd.pwdparms.get('Background',
                        [['chebyschev-1',True,3,1.0,0.0,0.0],{'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[],
                        'background PWDR':['',1.0,False]}]))
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Instrument Parameters'),
                    [Iparm1,Iparm2])
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Sample Parameters'),
                    rd.Sample)
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Peak List')
                    ,{'peaks':[],'sigDict':{}})
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Index Peak List'),
                    [[],[]])
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Unit Cells List'),
                    [])
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='Reflection Lists'),
                    {})
                # if any Control values have been set, move them into tree
                Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
                Controls.update(rd.Controls)
                # Tree entries complete
        # select and show last PWDR file to be read
        if Id:
            G2frame.EnablePlot = True
            G2frame.GPXtree.Expand(Id)
            G2frame.GPXtree.SelectItem(Id)
        dlg.Raise()
            
    def RunTimerGR(event):
        if GSASIIpath.GetConfigValue('debug'):
            import datetime
            print ("DBG_Timer tick at {:%d %b %Y %H:%M:%S}\n".format(datetime.datetime.now()))
        filelist = glob.glob(os.path.join(Settings['indir'],Settings['filter']))
        if not filelist: return
        #if GSASIIpath.GetConfigValue('debug'): print(filelist)
        Id = None
        for f in filelist:
            if f in Settings['filesread']: continue
            Settings['filesread'].append(f)
            rd = fileReaders[Settings['fmt']]
            rd.ReInitialize()
            if not rd.ContentsValidator(f):
                Settings['timer'].Stop()
                btnstart.SetLabel('Continue')
                G2MessageBox(dlg,'Error in reading file {}: {}'.format(
                    f, rd.errors))
                return
            #if len(rd.selections) > 1:
            #    G2fil.G2Print('Warning: Skipping file {}: multibank not yet implemented'.format(f))
            #    continue
            block = 0
            rdbuffer = {}
            repeat = True
            while repeat:
                repeat = False
                try:
                    flag = rd.Reader(f,buffer=rdbuffer, blocknum=block)
                except:
                    flag = False
                if flag:
                    rd.readfilename = f
                    if rd.warnings:
                        G2fil.G2Print("Read warning by", rd.formatName,
                                          "reader:",
                                          rd.warnings)
                    elif not block:
                        G2fil.G2Print("{} read by Reader {}"
                              .format(f,rd.formatName))
                    else:
                        G2fil.G2Print("{} block # {} read by Reader {}"
                                .format(f,block,rd.formatName))
                    block += 1    
                    repeat = rd.repeat
                else:
                    G2fil.G2Print("Warning: {} Reader failed to read {}"
                                      .format(rd.formatName,f))
                if 'phoenix' in wx.version():
                    HistName = 'PDF  '+rd.idstring
                else:
                    HistName = 'PDF  '+G2obj.StripUnicode(rd.idstring,'_')
                HistName = G2obj.MakeUniqueLabel(HistName,Settings['ReadList'])
                Settings['ReadList'].append(HistName)
                # put into tree
                Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=HistName)
                Ymin = np.min(rd.pdfdata[1])                 
                Ymax = np.max(rd.pdfdata[1])                 
                valuesdict = {
                    'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxsize),
                    'Offset':[0.0,0.0],'delOffset':0.02*Ymax,
                    'Yminmax':[Ymin,Ymax],
                    }
                G2frame.GPXtree.SetItemPyData(
                    G2frame.GPXtree.AppendItem(Id,text='PDF Controls'),
                    {'G(R)':[valuesdict,rd.pdfdata,HistName],
                         'diffGRname':'','diffMult':1.0,'Rmax':Ymax,})
                G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='PDF Peaks'),
                    {'Limits':[1.,5.],'Background':[2,[0.,-0.2*np.pi],False],'Peaks':[]})
                    
        # select and show last PWDR file to be read
        if Id:
            G2frame.EnablePlot = True
            G2frame.GPXtree.Expand(Id)
            G2frame.GPXtree.SelectItem(Id)
                
    global AutoLoadWindow    
    Settings = {}
    if AutoLoadWindow: # make sure only one window is open at a time
        try:
            AutoLoadWindow.Destroy()
        except:
            pass
        AutoLoadWindow = None
    if FileTyp == 'pwd':
        fileReaders = [i for i in G2fil.LoadImportRoutines("pwd", "Powder_Data")
                       if i.scriptable]
        fmtchoices = [p.longFormatName for p in fileReaders]
        Settings['fmt'] = [i for i,v in enumerate(fmtchoices) if 'fxye' in v][0]
    else:
        fileReaders = [i for i in G2frame.ImportPDFReaderlist]
#                       if i.scriptable]
        fmtchoices = [p.longFormatName for p in fileReaders]
        Settings['fmt'] = 0       
    Settings['ext'] = 0
    Settings['extStr'] = ''
    Settings['filter'] = '*.*'
    Settings['indir'] = os.getcwd()
    Settings['instfile'] = ''
    Settings['timer'] = wx.Timer()
    if FileTyp == 'pwd':
        Settings['timer'].Bind(wx.EVT_TIMER,RunTimerPWDR)
    else:
        Settings['timer'].Bind(wx.EVT_TIMER,RunTimerGR)
    Settings['filesread'] = []
    dlg = wx.Frame(G2frame,title='Automatic Data Loading',
                       style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
    mnpnl = wx.Panel(dlg)
    mnsizer = wx.BoxSizer(wx.VERTICAL)
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Select format:'))
    fmtSel = G2ChoiceButton(mnpnl,fmtchoices,Settings,'fmt',onChoice=onSetFmtSelection)
    sizer.Add(fmtSel)
    sizer.Add((-1,-1),1,wx.EXPAND,1)
    msg = '''This window serves two purposes. It can be used to read files 
as they are added to a directory or it can be used to read files from an 
externally-created file list. For either, set the file format and an 
instrument parameter file must be specified.
%%
* For automatic reading, the files must be found in the directory specified by
"Read from:" and the selected extension. The "File filter:" can be used to 
limit the files to those matching a wildcard, (for example, if 
"202408*pow*.*" is used as a filter, then files must begin with "202408" 
and must also contain the string "pow".) 
%%
* For reading from a list of files, press the "Read from file with a list 
of files" button. The input file must contain a list of files, one per line. 
Lines beginning in '#' are ignored. If more than one column is used 
(separated by commas or tabs), the file name should be the first column. 
File names can be in quotes, but this is not required. The extension 
is ignored, as is the "File filter". The "Read from" directory will be used 
if the file name does not contain a full path and the file is not in the 
current working directory.
'''
    sizer.Add(HelpButton(mnpnl,msg,wrap=400),0,wx.RIGHT,5)
    mnsizer.Add(sizer,0,wx.EXPAND)

    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Select extension:'))
    extSel = G2ChoiceButton(mnpnl,[],Settings,'ext',Settings,'extStr',onChoice=onSetExtSelection)
    sizer.Add(extSel,0)
    sizer.Add((-1,-1),1,wx.EXPAND,1)
    sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'  File filter: '))
    flterInp = ValidatedTxtCtrl(mnpnl,Settings,'filter')
    sizer.Add(flterInp)
    mnsizer.Add(sizer,0,wx.EXPAND,0)

    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Read from: '),0,wx.ALIGN_CENTER_VERTICAL)
    fInp3 = ValidatedTxtCtrl(mnpnl,Settings,'indir',size=(300,-1),OnLeave=TestInput)
    sizer.Add(fInp3,1,wx.EXPAND)
    btn3 = wx.Button(mnpnl,  wx.ID_ANY, "Browse")
    btn3.Bind(wx.EVT_BUTTON, OnBrowse)
    sizer.Add(btn3,0,wx.ALIGN_CENTER_VERTICAL)
    mnsizer.Add(sizer,0,wx.EXPAND)
        
    if FileTyp == 'pwd':
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(wx.StaticText(mnpnl, wx.ID_ANY,'Instrument parameter file from: '),0,wx.ALIGN_CENTER_VERTICAL)
        fInp4 = ValidatedTxtCtrl(mnpnl,Settings,'instfile',size=(300,-1),OnLeave=TestInput)
        sizer.Add(fInp4,1,wx.EXPAND)
        btn4 = wx.Button(mnpnl,  wx.ID_ANY, "Browse")
        btn4.Bind(wx.EVT_BUTTON, OnBrowse)
        sizer.Add(btn4,0,wx.ALIGN_CENTER_VERTICAL)
        mnsizer.Add(sizer,0,wx.EXPAND)
        # read a list of files
        FofFbtn = wx.Button(mnpnl,  wx.ID_ANY, 'Read from file with a list of files')
        FofFbtn.Bind(wx.EVT_BUTTON, OnFileOfFiles)
        mnsizer.Add(FofFbtn)
    # buttons on bottom
    sizer = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add((-1,-1),1,wx.EXPAND)
    btnstart = wx.Button(mnpnl,  wx.ID_ANY, "Start")
    btnstart.Bind(wx.EVT_BUTTON, OnStart)
    sizer.Add(btnstart)
    sizer.Add((20,-1),0,wx.EXPAND)
    btnclose = wx.Button(mnpnl,  wx.ID_ANY, "Close")
    onSetFmtSelection()
    btnclose.Bind(wx.EVT_BUTTON, OnQuit)
    sizer.Add(btnclose)
    sizer.Add((-1,-1),1,wx.EXPAND)
    mnsizer.Add(sizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP,5)
    mnpnl.SetSizer(mnsizer)
    mnsizer.Fit(dlg)
    dlg.CenterOnParent()
    dlg.Show()
    AutoLoadWindow = dlg # save window reference

# Deal with Origin 1/2 ambiguities ################################################################################
def ChooseOrigin(G2frame,rd):    
    G2elem.SetupGeneral(rd.Phase,G2frame.dirname)
    # make copy of Phase but shift atoms Origin 1->2
    O2Phase = copy.deepcopy(rd.Phase)
    # make copy of atoms, shift to alternate origin
    T = G2spc.spg2origins[rd.Phase['General']['SGData']['SpGrp']]
    O2atoms = O2Phase['Atoms']
    cx,ct,cs,cia = rd.Phase['General']['AtomPtrs']
    SGData = rd.Phase['General']['SGData']
    for atom in O2atoms:
        for i in [0,1,2]:
            atom[cx+i] += T[i]
            atom[cs:cs+2] = G2spc.SytSym(atom[3:6],SGData)[0:2] # update symmetry & mult
    #get density & distances
    DisAglData = {}
    DisAglData['SGData'] = rd.Phase['General']['SGData']
    DisAglData['Cell'] = rd.Phase['General']['Cell'][1:] #+ volume
    DisAglCtls = {'Factors': [0.85, 0],
                 'BondRadii': [], 'AngleRadii': [], 'AtomTypes': []}
    for atom in rd.Phase['Atoms']:
        DisAglCtls['BondRadii'].append(1.5)
        DisAglCtls['AngleRadii'].append(0)
        DisAglCtls['AtomTypes'].append(atom[ct])
    txt = ''
    for i,phObj in enumerate([rd.Phase,O2Phase]):
        if i:
            txt += "\n\nWith origin shift applied\n"
        else:
            txt += "\nWith current coordinates and original origin\n"
        cellContents = {}
        for atom in phObj['Atoms']:
            if atom[ct] in cellContents:
                cellContents[atom[ct]] += atom[cs+1]
            else:
                cellContents[atom[ct]] = atom[cs+1]
        txt += '   Unit cell Contents: '
        for i,k in enumerate(cellContents):
            if i: txt += ', '
            txt += '{}*{}'.format(cellContents[k],k)
        den,_ = G2mth.getDensity(phObj['General'])
        txt += "\n   Density {:.2f} g/cc\n".format(den)
                    
        DisAglData['OrigAtoms'] = DisAglData['TargAtoms'] = [
                        [i,]+atom[ct-1:ct+1]+atom[cx:cx+3] for
                        i,atom in enumerate(phObj['Atoms'])]
        # lbl,dis,angle = G2stMn.RetDistAngle(DisAglCtls,DisAglData)
        # # get unique distances
        # minDis = {} 
        # for i in dis: 
        #     for j,o,s,d,e in dis[i]: 
        #         key = '-'.join(sorted([lbl[i],lbl[j]])) 
        #         if key not in minDis: 
        #             minDis[key] = d 
        #         elif d < minDis[key]: 
        #             minDis[key] = d
        # thirdShortest = sorted([minDis[k] for k in minDis])[:3][-1]
        # shortTxt = ''
        # for k in minDis:
        #     if minDis[k] <= thirdShortest:
        #         if shortTxt: shortTxt += ', '
        #         shortTxt += "{}: {:.2f}".format(k,minDis[k])
        # txt += "   Shortest distances are "+shortTxt

    # do we know if there is a center of symmetry at origin?
    centro = None
    if 'xyz' in rd.SymOps:
        centro = False
        if '-x,-y,-z' in [i.replace(' ','').lower() for i in rd.SymOps['xyz']]:
            centro = True
            
    msg = 'Be careful here. This space group has two origin settings. GSAS-II requires the origin to be placed at a center of symmetry (Origin 2). You must choose the correct option below or all subsequent results will be *wrong*. For more info, press the help button (bottom right).\n'
    if centro:
        msg += '\nThere is an -x,-y,-z symmetry op in the file input, so this is likely already in Origin 2.\n'
    elif centro is None:
        msg += '\nNo symmetry operations are provided in the input file; you must review this yourself. You are recommended to review a plot of the structure to make sure the symmetry is correct.\n'
    else:
        msg += '\nSymmetry operations in the input file do not contain -x,-y,-z, indicating an origin shift is likely needed.\n'

    msg += '\nNote that the stoichometry computations below, made from the coordinates, may help indicate the correct origin choice:'

    width = 600
    dlg = wx.Dialog(G2frame,wx.ID_ANY,'Warning: Shift origin?',
                pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE,
                size=(width,-1))
    dlg.CenterOnParent()
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    txtbox = wx.StaticText(dlg,wx.ID_ANY,msg)
    txtbox.Wrap(width-10)
    mainSizer.Add(txtbox,0)
    mainSizer.Add((5,5))
    txtbox = wx.StaticText(dlg,wx.ID_ANY,txt)
    mainSizer.Add(txtbox,0,wx.ALIGN_CENTER,1)
    mainSizer.Add((10,10))

    O1Btn = wx.Button(dlg,wx.ID_ANY,"Apply origin shift")
    O1Btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_OK))
    O2Btn = wx.Button(dlg,wx.ID_ANY,"Keep current coordinates")
    O2Btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_YES))
    if centro:
        O2Btn.SetDefault()
    elif centro is not None:
        O1Btn.SetDefault()
    btnSizer = wx.BoxSizer(wx.HORIZONTAL)
    btnSizer.Add((20,20),1)
    btnSizer.Add(O1Btn)
    btnSizer.Add((10,10),0)
    btnSizer.Add(O2Btn)
    btnSizer.Add((20,20),1)
    btnSizer.Add(HelpButton(dlg,helpIndex='Origin1'),0,wx.RIGHT,5)
    mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
    dlg.SetSizer(mainSizer)
    dlg.Fit()
    ans = dlg.ShowModal()
    if ans == wx.ID_OK:
        dlg.Destroy()
        return O2Phase
    elif ans == wx.ID_YES:
        dlg.Destroy()
        return rd.Phase
    else:
        dlg.Destroy()
        return None

def makeContourSliders(G2frame,Ymax,PlotPatterns,newPlot,plottype):
    '''Create a non-modal dialog for sliders to set contour plot 
    intensity thresholds. 
    '''
    def updatePlot():
        'updates plot after a change in values'
        wx.CallAfter(PlotPatterns,G2frame,newPlot=newPlot,plotType=plottype)
    def OnSlider(event):
        'respond when min or max slider is moved'
        obj = event.GetEventObject()
        val = obj.GetValue()/100.
        if obj.mode == 'max':
            G2frame.Cmax = val
        else:
            G2frame.Cmin = val
        obj.txt.SetValue(int(Ymax*val))
        updatePlot()
    def OnNewVal(*args,**kwargs):
        'respond when a value is placed in the min or max text box'
        obj = kwargs['tc']
        if obj.mode == 'max':
            val = Range[1]
            G2frame.Cmax = val/Ymax
        else:
            val = Range[0]
            G2frame.Cmin = val/Ymax
        obj.slider.SetValue(int((100*val/Ymax) + 0.5))
        updatePlot()
    # makeContourSliders starts here
    Range = [Ymax*G2frame.Cmin,Ymax*G2frame.Cmax]
    dlg = wx.Dialog(G2frame.plotFrame,style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
    vbox = wx.BoxSizer(wx.VERTICAL)
    vbox.Add((-1,5))
    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add((-1,-1),1,wx.EXPAND,1)
    hbox.Add(wx.StaticText(
                dlg,wx.ID_ANY,'Set Contour Intensity Limits'),1,wx.ALL)
    hbox.Add((-1,-1),1,wx.EXPAND,1)
    vbox.Add(hbox)
    vbox.Add((-1,10))
    dlg.slideSizer = wx.FlexGridSizer(2,3,5,5)
    dlg.slideSizer.AddGrowableCol(2)

    dlg.slideSizer.Add(wx.StaticText(parent=dlg,label=' Min intensity'),0,WACV)
    minSel = wx.Slider(parent=dlg,style=wx.SL_HORIZONTAL,
                           value=int(100*G2frame.Cmin+0.5))
    minSel.Bind(wx.EVT_SLIDER, OnSlider)
    minVal = ValidatedTxtCtrl(dlg,Range,0,xmin=0,
                xmax=Ymax-1, OnLeave=OnNewVal)
    minVal.slider = minSel
    minVal.mode = 'min'
    minSel.txt = minVal
    minSel.mode = 'min'
    dlg.slideSizer.Add(minVal,0,WACV)
    dlg.slideSizer.Add(minSel,0,wx.EXPAND|wx.ALL,3)

    dlg.slideSizer.Add(wx.StaticText(parent=dlg,label=' Max intensity'),0,WACV)
    maxSel = wx.Slider(parent=dlg,style=wx.SL_HORIZONTAL,
                           value=int(100*G2frame.Cmax+0.5))
    maxSel.Bind(wx.EVT_SLIDER, OnSlider)
    maxVal = ValidatedTxtCtrl(dlg,Range,1,xmin=1,
                xmax=Ymax, OnLeave=OnNewVal)
    maxVal.slider = maxSel
    maxVal.mode = 'max'
    maxSel.txt = maxVal
    maxSel.mode = 'max'
    dlg.slideSizer.Add(maxVal,0,WACV)
    dlg.slideSizer.Add(maxSel,0,wx.EXPAND|wx.ALL,3)

    vbox.Add(dlg.slideSizer,0,wx.EXPAND,0)
    vbox.Add((-1,15))

    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add((-1,-1),1,wx.EXPAND,1)
    btn = wx.Button(dlg, wx.ID_CLOSE) 
    hbox.Add(btn,0,wx.ALL,0)
    btn.Bind(wx.EVT_BUTTON,lambda x: dlg.Destroy())
    hbox.Add((-1,-1),1,wx.EXPAND,1)
    vbox.Add(hbox,0,wx.EXPAND,0)
    vbox.Add((-1,5))

    dlg.SetSizer(vbox)
    vbox.Fit(dlg)
    wx.CallLater(100,minVal.SetValue,Range[0])
    dlg.Show()

################################################################################
# GPX browser routines
def skimGPX(fl):
    '''pull out fit information from a .gpx file quickly

    :returns: dict with status info
    '''
    if fl is None: return {}
    result = {'other':[]}
    if not os.path.exists(fl):
        return {'error':'File does not exist!'}
    cnt = 0
    hist = 0
    fp = open(fl,'rb')
    result['last saved'] = time.ctime(os.stat(fl).st_mtime)
    try:
        while True:
            cnt += 1
            note = None
            try:
                data = G2IO.cPickleLoad(fp)
            except EOFError:
                #print(cnt,'entries read')        
                break
            if cnt > 50:  # don't spend too long on this file, if big
                result['PWDR'] += 3*['   .']
                break
            datum = data[0]
            if datum[0] == 'Notebook':
                result[datum[0]] = datum[1][-1]
            elif 'Controls' in datum[0]:
#                datum[0]['Seq Data']
                if 'LastSavedUsing' in datum[1]:
                    result['last saved'] += ' (v' + datum[1]['LastSavedUsing'] +')'
                    
            elif datum[0] == 'Covariance':
                d = datum[1].get('Rvals')
                if d:
                    result[datum[0]] = 'Overall: Rwp={:.2f}, GOF={:.1f}'.format(
                        d.get('Rwp','?'),d.get('GOF','?'))
                    if d.get('converged',False): result[datum[0]] += '  **Converged**'
            elif datum[0].startswith('PWDR '):
                if 'Residuals' not in datum[1][0]: continue
                if 'PWDR' not in result: result['PWDR'] = []
                result['PWDR'].append(
                    "hist #{}: wR={:.2f} ({:})".format(
                        hist,datum[1][0]['Residuals'].get('wR','?'),datum[0]))
                hist += 1
            elif datum[0].startswith('HKLF '):
                note = 'Single crystal histogram(s)'
            elif datum[0].startswith('REFD '):
                note = 'Reflectivity histogram(s)'
            elif datum[0].startswith('SASD '):
                note = 'Small angle histogram(s)'
            elif datum[0].startswith('PDF  '):
                note = 'PDF histogram(s)'
            elif datum[0].startswith('IMG '):
                note = 'Image(s)'
            elif datum[0] == 'Sequential results':
                note = 'Sequential results'
            elif datum[0] in ('Constraints','Restraints','Rigid bodies'):
                pass
            else:
#                print(datum[0])
#                breakpoint()
                pass   
            if note:
                if note not in result['other']:
                    result['other'].append(note)
    except Exception as msg:
        result['error'] = 'read error: '+str(msg)
    finally:
        fp.close()
    return result

class gpxFileSelector(wx.Dialog):
    '''Create a file selection widget for locating .gpx files as a modal
    dialog. Displays status information on selected files. After creating 
    this use dlg.ShowModal() to wait for selection of a file.
    If dlg.ShowModal() returns wx.ID_OK, use dlg.Selection (multiple=False)
    to obtain the selected file or dlg.Selections (multiple=True) to 
    obtain a list of multiple files.

    :param wx.Frame parent: name of panel or frame that will be
      the parent to the dialog. Can be None.

    :param path startdir: Specifies the initial directory that is
      opened when the window is initially opened. Default is '.'

    :param bool multiple: if True, checkboxes are used to allow
      selection of multiple files. Default is False

    '''
    def __init__(self,parent,startdir='.',multiple=False,*args,**kwargs):
        self.timer = None
        self.delay = 1500 # time to wait before applying filter (1.5 sec)
        self.Selection = None
        self.Selections = []
        self.startDir = startdir
        if startdir == '.':
            self.startDir = os.getcwd()
        self.multiple = multiple
        wx.Dialog.__init__(self, parent=parent, 
                                 style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.CenterOnParent()
        
        topSizer = wx.BoxSizer(wx.VERTICAL)
        self.dirBtn = wxfilebrowse.DirBrowseButton(self,wx.ID_ANY, size=(650, -1), 
                            changeCallback = self.DirSelected,
                            startDirectory = self.startDir
                    )
        topSizer.Add(self.dirBtn,0,wx.EXPAND,1)
        
        subSiz = wx.BoxSizer(wx.HORIZONTAL)
        self.opt = {'useBak':False, 'sort':0, 'filter':'*'}
        chk = G2CheckBoxFrontLbl(self,' Include .bakXX?',self.opt,'useBak',
                                     OnChange=self.DirSelected)
        subSiz.Add(chk,0,wx.ALIGN_CENTER_VERTICAL,0)
        subSiz.Add((10,-1),1,wx.EXPAND,1)
        subSiz.Add(wx.StaticText(self,wx.ID_ANY,'   Sort by: '),0,wx.ALIGN_CENTER_VERTICAL,1)
        choices = ['age','name (alpha+case)','name (alpha)']
        for w in G2RadioButtons(self,self.opt,'sort',choices,
                                    OnChange=self.DirSelected):
            subSiz.Add(w,0,wx.ALIGN_CENTER_VERTICAL,0)
        subSiz.Add((10,-1),1,wx.EXPAND,1)
        subSiz.Add(wx.StaticText(self,wx.ID_ANY,'Name \nFilter: '),0,wx.ALIGN_CENTER_VERTICAL,1)
        self.filterBox = ValidatedTxtCtrl(self, self.opt, 'filter', 
                                size=(80,-1), style=wx.TE_PROCESS_ENTER,
                                OnLeave=self.DirSelected, notBlank=False)
        self.filterBox.Bind(wx.EVT_TEXT,self._startUpdateTimer)
        self.filterBox.Bind(wx.EVT_TEXT_ENTER,self.DirSelected)
        subSiz.Add(self.filterBox)
        subSiz.Add((2,-1))
        
        topSizer.Add(subSiz,0,wx.EXPAND,0)

        mainPanel = wx.SplitterWindow(self, wx.ID_ANY, style=wx.SP_LIVE_UPDATE|wx.SP_3D)
        mainPanel.SetMinimumPaneSize(100)

        if self.multiple:
            self.fileBox = wx.CheckListBox(mainPanel,wx.ID_ANY, size=(200, 200),
                                               style=wx.LB_SINGLE)
            self.fileBox.Bind(wx.EVT_CHECKLISTBOX,self.FileSelected)
        else:
            self.fileBox = wx.ListBox(mainPanel,wx.ID_ANY, size=(200, 200),
                                          style=wx.LB_SINGLE)
        self.fileBox.Bind(wx.EVT_LISTBOX,self.FileSelected)

        self.rtc = wxrt.RichTextCtrl(mainPanel, style=wx.VSCROLL|wx.HSCROLL|
                                       wx.NO_BORDER|wx.richtext.RE_READONLY)
        mainPanel.SplitVertically(self.fileBox, self.rtc, 200)
        topSizer.Add(mainPanel,1,wx.EXPAND)
        
        subSiz = wx.BoxSizer(wx.HORIZONTAL)
        subSiz.Add((-1,-1),1,wx.EXPAND,1)
        self.OKbtn = wx.Button(self, wx.ID_OK, label='Open')
        self.OKbtn.Enable(False)   # A file must be selected 1st
        btn = wx.Button(self, wx.ID_CANCEL)
        subSiz.Add(self.OKbtn)
        subSiz.Add((5,-1))
        subSiz.Add(btn)
        subSiz.Add((-1,-1),1,wx.EXPAND,1)
        topSizer.Add((-1,5))
        topSizer.Add(subSiz,0,wx.EXPAND)
        topSizer.Add((-1,5))
        self.SetSizer(topSizer)
        topSizer.Fit(self)
        self.dirBtn.SetValue(self.startDir)

    def _startUpdateTimer(self,event):
        if self.timer:
            self.timer.Restart(self.delay)
        else:
            self.timer = wx.CallLater(self.delay,self.DirSelected)
        
    def DirSelected(self,event=None,*args,**kwargs):
        '''Respond to a directory being selected. List files found in fileBox and 
        clear any selections. Also clear any reference to a timer. 
        '''
        import re
        try:
            if self.timer: self.timer.Stop()
        except:
            pass
        self.timer = None
        self.fileBox.Clear()
        self.rtc.Clear()
        self.Selection = None
        self.Selections = []
        self.OKbtn.Enable(False)
        glb = self.opt['filter'].strip()
        if not glb:
            glb = '*'
        elif not '*' in glb:
            glb = '*' + glb + '*'
        fullglob = os.path.join(self.dirBtn.GetValue(),glb+'.gpx')
        self.fl = glob.glob(fullglob)
        if self.opt['useBak']:
            self.sl = [(os.path.split(i)[1],os.stat(i).st_mtime,i) for i in self.fl]
        else:
            self.sl = [(os.path.split(i)[1],os.stat(i).st_mtime,i) for i in self.fl 
                  if not re.match(r'.*\.bak\d+\.gpx.*',i)]
        if self.opt['sort'] == 0:
            self.sl.sort(key=lambda x: x[1],reverse=True)
        elif self.opt['sort'] == 1:
            self.sl.sort(key=lambda x: x[0])
        else:
            self.sl.sort(key=lambda x: x[0].lower())
        items = [i[0]+' ('+self._fmtTimeStampDelta(i[1])+')' for i in self.sl]
        if items: 
            self.fileBox.InsertItems(items,0)
        
    def FileSelected(self,event): 
        '''Respond to a file being selected (or checked in multiple mode)
        '''
        if self.multiple:  # disable  Open when nothing is selected
            self.Selections = []
            OK = False
            for i in self.fileBox.GetCheckedItems():
                self.Selections.append(self.sl[i][2])
                OK = True
            self.OKbtn.Enable(OK)
        else:
            self.OKbtn.Enable(True)
        self.Selection = self.sl[self.fileBox.GetSelection()][2]
        result = skimGPX(self.Selection)
        self.displayGPXrtc(result,self.Selection)

    def displayGPXrtc(self,result,fwp):
        '''Show info about selected file in a RichText display'''
        self.rtc.Clear()
        if fwp is None: return
        self.rtc.Freeze()
        self.rtc.BeginSuppressUndo()
        self.rtc.BeginAlignment(wx.TEXT_ALIGNMENT_CENTER)
        self.rtc.BeginFontSize(14)
        self.rtc.BeginBold()
        self.rtc.WriteText(os.path.split(fwp)[1])
        self.rtc.EndBold()
        self.rtc.Newline()
        self.rtc.EndFontSize()
        self.rtc.EndAlignment()
        self.rtc.WriteText('last saved on ')
        self.rtc.WriteText(result['last saved'])
        self.rtc.Newline()
        if 'Covariance' in result:
            self.rtc.BeginLeftIndent(0,40)
            self.rtc.WriteText(result['Covariance'])
            self.rtc.Newline()
            self.rtc.EndLeftIndent()
        if 'Notebook' in result and len(result.get('Notebook','').strip()):
            self.rtc.BeginLeftIndent(0,40)
            self.rtc.BeginItalic()
            self.rtc.WriteText('Last notebook entry: ')
            self.rtc.EndItalic()
            self.rtc.WriteText(result['Notebook'])
            self.rtc.Newline()
            self.rtc.EndLeftIndent()
        if len(result.get('other',[])) > 0:
            self.rtc.BeginParagraphSpacing(0,0)
            self.rtc.BeginLeftIndent(0)
            self.rtc.BeginBold()
            self.rtc.WriteText('Data types in project:')
            self.rtc.EndBold()
            self.rtc.EndLeftIndent()
            self.rtc.Newline()
            self.rtc.BeginLeftIndent(40)
            for line in result['other']:
                self.rtc.WriteText(line+'\n')
            self.rtc.EndLeftIndent()
            self.rtc.EndParagraphSpacing()
                   
        if 'PWDR' in result:
            self.rtc.BeginParagraphSpacing(0,0)
            self.rtc.BeginLeftIndent(0)
            self.rtc.BeginBold()
            self.rtc.WriteText('Powder histograms:')
            self.rtc.EndBold()
            self.rtc.EndLeftIndent()
            self.rtc.Newline()
            self.rtc.BeginLeftIndent(40)
            for line in result['PWDR']:
                self.rtc.WriteText(line+'\n')
            self.rtc.EndLeftIndent()
            self.rtc.EndParagraphSpacing()
            
        if 'error' in result:
            self.rtc.Newline()
            self.rtc.BeginBold()
            self.rtc.WriteText('Error encountered: ')
            self.rtc.EndBold()
            self.rtc.WriteText(result['error'])
        self.rtc.EndSuppressUndo()
        self.rtc.Thaw()

    def _fmtTimeStampDelta(self,tm):
        'Show file age relative to now'
        delta = time.time() - tm
        if delta > 60*60*24*365:
            return "{:.2f} years".format(delta/(60*60*24*365))
        elif delta > 60*60*24*7:
            return "{:.1f} weeks".format(delta/(60*60*24*7))
        elif delta > 60*60*24:
            return "{:.1f} days".format(delta/(60*60*24))
        elif delta > 60*60:
            return "{:.1f} hours".format(delta/(60*60))
        else:
            return "{:.1f} minutes".format(delta/60)

def setColorButton(parent,array,key,callback=None,callbackArgs=[]):
    '''Define a button for setting colors
    This bypasses the bug in wx4.1.x in ColourSelect
    '''
    import wx.lib.colourselect as wcs
    def OnColor(event):
        array[key] = list(event.GetValue())[:3]
        if callback: callback(*callbackArgs)
    def onSetColour(event):
        dlg = wx.ColourDialog(parent.GetTopLevelParent())
        try:
            dlg.GetColourData().SetChooseFull(True)
            dlg.GetColourData().SetColour(array[key])
            if dlg.ShowModal() == wx.ID_OK:
                data = dlg.GetColourData()
                array[key] = data.GetColour().Get()
                print('OK',array[key])
        finally:
            dlg.Destroy()    
    if wx.__version__.startswith('4.1'):
        colorButton = wx.Button(parent,wx.ID_ANY,'Set')
        colorButton.Bind(wx.EVT_BUTTON, onSetColour)
    else:
        colorButton = wcs.ColourSelect(parent,colour=array[key],size=wx.Size(25,25))
        colorButton.Bind(wcs.EVT_COLOURSELECT, OnColor)
    return colorButton

def NISTlatUse(msgonly=False):
        msg = '''Performing cell symmetry search using NIST*LATTICE. Please cite:
        V. L. Karen and A. D. Mighell, NIST Technical Note 1290 (1991),
        https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1290.pdf  
        and 
        V. L. Karen & A. D. Mighell, U.S. Patent 5,235,523,
        https://patents.google.com/patent/US5235523A/en?oq=5235523'''
        print(msg)
        if msgonly: return msg
        wx.MessageBox(msg,caption='Using NIST*LATTICE',style=wx.ICON_INFORMATION)

def Load2Cells(G2frame,phase):
    '''Accept two unit cells and use NIST*LATTICE to search for a relationship 
    that relates them. 

    The first unit cell is initialized as the currently selected phase and
    the second unit cell is set to the first different phase from the tree. 
    The user can initialize the cell parameters to select a different phase 
    for either cell or can type in the values themselves. 

    :param wx.Frame G2frame: The main GSAS-II window
    :param dict phase: the currently selected frame
    '''
    def setRatioMax(*arg,**kwarg):
        '''Set the value for the max volume used in the search according 
        to the type of search selected. 
        '''
        if nistInput[2] == 'I':
            volRatW.Validator.xmax = 40
            volMaxLbl.SetLabel(' (ratio, 1 to 40)')
        else:
            volRatW.Validator.xmax = 10
            volMaxLbl.SetLabel(' (ratio, 1 to 10)')
        volRatW.SetValue(min(volRatW.Validator.xmax,nistInput[3]))
    def computeNISTlatCompare(event):
        'run NIST*LATTICE after the compute button is pressed'
        import nistlat
        out = nistlat.CompareCell(cellLen[0], cellCntr[0],
                                  cellLen[1], cellCntr[1],
                    tolerance=3*[nistInput[0]]+3*[nistInput[1]],
                    mode=nistInput[2], vrange=nistInput[3])
        if len(out):
            msg = str(len(out))+' Transformations were found. See console for matrices.'
            G2MessageBox(G2frame,msg,'Transforms found')
            print(len(out),'transform matrices found')
            for i in out:
                print(' ',i[5][0],'/',i[5][1],'/',i[5][2])
        else:
            G2MessageBox(G2frame,
                'No transforms were found within supplied limits',
                'No transforms found')
    def setCellFromPhase(event):
        '''respond to "set from phase" button. A phase is selected and 
        the unit cell info is loaded from that phase into the appropriate
        cell widgets.
        '''
        cell = event.GetEventObject().cellNum
        widgets = event.GetEventObject().widgets
        phaseList = list(Phases.keys())
        if len(Phases) == 0:
            print('No phases in project')
            return
        elif len(Phases) == 1:
            p = phaseList[0]
        else:
            dlg = G2SingleChoiceDialog(G2frame,'Select a phase from list',
                                           'Select phase',phaseList)
            dlg.CenterOnParent()
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    p = phaseList[dlg.GetSelection()]
                else:
                    return
            finally:
                dlg.Destroy()
        ph2 = Phases[p]
        cellLen[cell] = ph2['General']['Cell'][1:7]
        cellCntr[cell] = ph2['General']['SGData']['SpGrp'].strip()[0]
        for val,wid in zip(cellLen[cell],widgets[:6]):
            wid.SetValue(val)
        widgets[6].SetValue(cellCntr[cell])
        dlg.Raise() # needed to bring modal dialog to front, at least on Mac
        
    # Load2Cells starts here  
    msg = NISTlatUse(True)
    nistInput=[0.2,1.,'I',8]
    cellLen = [None,None]
    cellCntr = [None,None]
    cellLen[0] = phase['General']['Cell'][1:7]
    cellCntr[0] = phase['General']['SGData']['SpGrp'].strip()[0]
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    cellLen[1] = 3*[1.]+3*[90.]
    cellCntr[1] = 'P'
    for p in Phases:
        ph2 = Phases[p]
        if ph2['ranId'] != phase['ranId']:
            cellLen[1] = ph2['General']['Cell'][1:7]
            cellCntr[1] = ph2['General']['SGData']['SpGrp'].strip()[0]
            break

    dlg = wx.Dialog(G2frame,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    dlg.CenterOnParent()

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(wx.StaticText(dlg,label='NIST*LATTICE: Relate Two Unit Cells'),
                  0,wx.ALIGN_CENTER_HORIZONTAL,0)
    sizer.Add((-1,15))
    sizer.Add(wx.StaticText(dlg,label=msg))
    sizer.Add((-1,15))
    tableSizer = wx.FlexGridSizer(0,9,0,0)
    tableSizer.Add((-1,-1))
    for l in u'abc\u03B1\u03B2\u03B3':
        tableSizer.Add(wx.StaticText(dlg,label=l),0,WACV|wx.ALIGN_CENTER)
    tableSizer.Add(wx.StaticText(dlg,label='Centering'),0,WACV|wx.ALIGN_LEFT)
    tableSizer.Add((-1,-1))
    for cell in range(2):
        tableSizer.Add(wx.StaticText(dlg,label='Cell '+str(cell+1)),0,wx.ALIGN_CENTER|wx.RIGHT,5)
        wlist = []
        for i in range(6):
            l = 3
            if i < 3: l = 4
            w = ValidatedTxtCtrl(dlg,cellLen[cell],i,nDig=(7,l),size=(60,-1))
            wlist.append(w)
            tableSizer.Add(w,0,wx.LEFT|wx.RIGHT,3)
        w = EnumSelector(dlg,cellCntr,cell,['P','A','B','C','F','I','R'])
        tableSizer.Add(w)
        wlist.append(w)
        btn = wx.Button(dlg, wx.ID_ANY, 'Load from phase')
        btn.cellNum = cell
        btn.widgets = wlist
        btn.Bind(wx.EVT_BUTTON, setCellFromPhase)
        tableSizer.Add(btn)
    sizer.Add(tableSizer,0,wx.LEFT|wx.RIGHT,20)
    tableSizer = wx.FlexGridSizer(0,3,0,0)
    tableSizer.Add(wx.StaticText(dlg,label='Cell length tolerance: '),
                       0,WACV|wx.ALIGN_LEFT)
    w = ValidatedTxtCtrl(dlg,nistInput,0,nDig=(6,2))
    tableSizer.Add(w)
    tableSizer.Add(wx.StaticText(dlg,label=' (A) '))
    tableSizer.Add(wx.StaticText(dlg,label='Cell angle tolerance: '),
                       0,WACV|wx.ALIGN_LEFT)
    w = ValidatedTxtCtrl(dlg,nistInput,1,nDig=(6,1))
    tableSizer.Add(w)
    tableSizer.Add(wx.StaticText(dlg,label=' (degrees) '))
    tableSizer.Add(wx.StaticText(dlg,label='Cell volume range: '),
            0,WACV|wx.ALIGN_LEFT)
    volRatW = ValidatedTxtCtrl(dlg,nistInput,3,xmin=1,xmax=40)
    tableSizer.Add(volRatW)
    volMaxLbl = wx.StaticText(dlg,label=' (ratio, 1 to 40)')
    tableSizer.Add(volMaxLbl)
    sizer.Add(tableSizer,0,wx.EXPAND)
    tableSizer = wx.FlexGridSizer(0,2,0,0)
    tableSizer.Add(wx.StaticText(dlg,label='Search mode: '),0,WACV|wx.ALIGN_LEFT)
    tableSizer.Add(EnumSelector(dlg,nistInput,2,['Integral matrices', 'Fractional matrices'],
        ['I','F'],OnChange=setRatioMax))
    sizer.Add(tableSizer,0,wx.EXPAND)
    btnSizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = wx.Button(dlg, wx.ID_ANY,'Compute')
    btn.Bind(wx.EVT_BUTTON, computeNISTlatCompare)
    btnSizer.Add((-1,-1),1,wx.EXPAND)
    btnSizer.Add(btn)
    btnSizer.Add((-1,-1),1,wx.EXPAND)
    sizer.Add(btnSizer,0,wx.EXPAND|wx.CENTER,0)

    btnsizer = wx.StdDialogButtonSizer()
    btn = wx.Button(dlg, wx.ID_CLOSE)
    btn.SetDefault()
    btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_OK))
    btnsizer.AddButton(btn)
    btnsizer.Realize()
    sizer.Add(btnsizer, 0, wx.EXPAND|wx.ALL, 5)
    dlg.SetSizer(sizer)
    sizer.Fit(dlg)
    
    if dlg.ShowModal() == wx.ID_OK:
        dlg.Destroy()    
    else:
        dlg.Destroy()    
        return

class ScrolledStaticText(wx.StaticText):
    '''Fits a long string into a small space by scrolling it. Inspired by 
    ActiveText.py from J Healey <rolfofsaxony@gmx.com> 
    https://discuss.wxpython.org/t/activetext-rather-than-statictext/36370

    Use examples::

      frm = wx.Frame(None) # create a frame
      ms = wx.BoxSizer(wx.VERTICAL)
      text = 'this is a long string that will be scrolled'
      ms.Add(G2G.ScrolledStaticText(frm,label=text))
      txt = G2G.ScrolledStaticText(frm,label=text, lbllen=20)
      smallfont = wx.SystemSettings.GetFont(wx.SYS_SYSTEM_FONT)
      smallfont.SetPointSize(10)
      txt.SetFont(smallfont)
      ms.Add(txt)
      ms.Add(G2G.ScrolledStaticText(frm,label=text,dots=False,delay=250,lbllen=20))
      frm.SetSizer(ms)
    
    :param w.Frame parent: Frame or Panel where widget will be placed
    :param str label: string to be displayed
    :param int delay: time between updates in ms (default is 100)
    :param int lbllen: number of characters to show (default is 15)
    :param bool dots: If True (default) ellipsis (...) are placed 
        at the beginning and end of the string when any characters 
        in the string are not shown. The displayed string length 
        will thus be lbllen+6 most of the time
    :param (other): other optional keyword parameters for the
      wx.StaticText widget such as size or style may be specified.
    '''
    def __init__(self, parent, label='', delay=100, lbllen=15, 
                 dots=True, **kwargs):
        wx.StaticText.__init__(self, parent, wx.ID_ANY, '', **kwargs)
        self.fullmsg = label
        self.lbllen = lbllen
        self.msgpos = 0
        self.dots = dots
        self.onTimer(None)
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.onTimer)
        self.timer.Start(delay, wx.TIMER_CONTINUOUS)

    def onTimer(self,event):
        if self.dots and self.msgpos > 0:
            txt = '...'
        else:
            txt = ''
        txt += self.fullmsg[self.msgpos:self.msgpos+self.lbllen+1]
        if self.dots and self.msgpos+self.lbllen < len(self.fullmsg):
            txt += '...'
        try:
            self.SetLabel(txt)
        except:
            self.timer.Stop()
        self.msgpos += 1
        if self.msgpos >= len(self.fullmsg): self.msgpos = 0

#===========================================================================
# this has been moved to GSASIIfiles, since it does not need wx
# def openInNewTerm(project=None,g2script=None,pythonapp=sys.executable):
#     '''Open a new and independent GSAS-II session in separate terminal 
#     or console window and as a separate process that will continue
#     even if the calling process exits.
#     Intended to work on all platforms. 

#     This could be used to run other scripts inside python other than GSAS-II

#     :param str project: the name of an optional parameter to be
#       passed to the script (usually a .gpx file to be opened in 
#       a new GSAS-II session)
#     :param str g2script: the script to be run. If None (default)
#       the GSASII.py file in the same directory as this file will
#       be used. 
#     :param str pythonapp: the Python interpreter to be used. 
#       Defaults to sys.executable which is usually what is wanted.
#     :param str terminal: a name for a preferred terminal emulator
#     '''
#     import subprocess
#     if g2script is None:
#         g2script = os.path.join(os.path.dirname(__file__),'GSASII.py')
    
#     if sys.platform == "darwin":
#         if project:
#             script = f'''
# set python to "{pythonapp}"
# set appwithpath to "{g2script}"
# set filename to "{project}"
# set filename to the quoted form of the POSIX path of filename

# tell application "Terminal"
#      activate
#      do script python & " " & appwithpath & " " & filename & "; exit"
# end tell
# '''
#         else:
#             script = f'''
# set python to "{pythonapp}"
# set appwithpath to "{g2script}"

# tell application "Terminal"
#      activate
#      do script python & " " & appwithpath & " " & "; exit"
# end tell
# '''
#         subprocess.Popen(["osascript","-e",script])
#     elif sys.platform.startswith("win"):
#         cmds = [pythonapp, g2script]
#         if project: cmds += [project]
#         subprocess.Popen(cmds,creationflags=subprocess.CREATE_NEW_CONSOLE)
#     else:
#         import shutil
#         script = ''
#         # try a bunch of common terminal emulators in Linux
#         # there does not appear to be a good way to way to specify this
#         # perhaps this should be a GSAS-II config option
#         for term in ("lxterminal", "gnome-terminal", 'konsole', "xterm",
#                          "terminator", "terminology", "tilix"):
#             try:
#                 found = shutil.which(term)
#                 if not found: continue
#             except AttributeError:
#                 print(f'shutil.which() failed (why?); assuming {term} present')
#                 found = True
#             if term == "gnome-terminal":
#                 #terminal = 'gnome-terminal -t "GSAS-II console" --'
#                 cmds = [term,'--title','"GSAS-II console"','--']
#                 script = "echo; echo Press Enter to close window; read line"
#                 break
#             elif term == "lxterminal":
#                #terminal = 'lxterminal -t "GSAS-II console" -e'
#                cmds = [term,'-t','"GSAS-II console"','-e']
#                script = "echo;echo Press Enter to close window; read line"
#                break
#             elif term == "xterm":
#                 #terminal = 'xterm -title "GSAS-II console" -hold -e'
#                 cmds = [term,'-title','"GSAS-II console"','-hold','-e']
#                 script = "echo; echo This window can now be closed"
#                 break
#             elif term == "terminator":
#                 cmds = [term,'-T','"GSAS-II console"','-x']
#                 script = "echo;echo Press Enter to close window; read line"
#                 break
#             elif term == "konsole":
#                 cmds = [term,'-p','tabtitle="GSAS-II console"','--hold','-e']
#                 script = "echo; echo This window can now be closed"
#                 break
#             elif term == "tilix":
#                 cmds = [term,'-t','"GSAS-II console"','-e']
#                 script = "echo;echo Press Enter to close window; read line"
#                 break
#             elif term == "terminology":
#                 cmds = [term,'-T="GSAS-II console"','--hold','-e']
#                 script = "echo; echo This window can now be closed"
#                 break                
#         else:
#             print("No known terminal was found to use, Can't start {}")
#             return

#         fil = '/tmp/GSAS2-launch.sh'
#         cmds += ['/bin/sh',fil]
#         fp = open(fil,'w')
#         if project:
#             fp.write(f"{pythonapp} {g2script} {project}\n")
#         else:
#             fp.write(f"{pythonapp} {g2script}\n")
#         fp.write(f"rm {fil}\n")
#         if script:
#             fp.write(f"{script}\n")
#         fp.close()
#         subprocess.Popen(cmds,start_new_session=True)

#===========================================================================
def ExtractFileFromZip(filename, selection=None, confirmread=True,
                       confirmoverwrite=True, parent=None,
                       multipleselect=False):
    '''If the filename is a zip file, extract a file from that
    archive.

    :param list Selection: used to predefine the name of the file
      to be extracted. Filename case and zip directory name are
      ignored in selection; the first matching file is used.

    :param bool confirmread: if True asks the user to confirm before expanding
      the only file in a zip

    :param bool confirmoverwrite: if True asks the user to confirm
      before overwriting if the extracted file already exists

    :param bool multipleselect: if True allows more than one zip
      file to be extracted, a list of file(s) is returned.
      If only one file is present, do not ask which one, otherwise
      offer a list of choices (unless selection is used).
    
    :returns: the name of the file that has been created or a
      list of files (see multipleselect)

    If the file is not a zipfile, return the name of the input file.
    If the zipfile is empty or no file has been selected, return None
    '''
    import zipfile # do this now, since we can save startup time by doing this only on need
    import shutil
    zloc = os.path.split(filename)[0]
    if not zipfile.is_zipfile(filename):
        #print("not zip")
        return filename

    z = zipfile.ZipFile(filename,'r')
    zinfo = z.infolist()

    if len(zinfo) == 0:
        #print('Zip has no files!')
        zlist = [-1]
    if selection:
        choices = [os.path.split(i.filename)[1].lower() for i in zinfo]
        if selection.lower() in choices:
            zlist = [choices.index(selection.lower())]
        else:
            print('debug: file '+str(selection)+' was not found in '+str(filename))
            zlist = [-1]
    elif len(zinfo) == 1 and confirmread:
        result = wx.ID_NO
        dlg = wx.MessageDialog(
            parent,
            'Is file '+str(zinfo[0].filename)+
            ' what you want to extract from '+
            str(os.path.split(filename)[1])+'?',
            'Confirm file', 
            wx.YES_NO | wx.ICON_QUESTION)
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result == wx.ID_NO:
            zlist = [-1]
        else:
            zlist = [0]
    elif len(zinfo) == 1:
        zlist = [0]
    elif multipleselect:
        # select one or more from a from list
        choices = [i.filename for i in zinfo]
        dlg = G2MultiChoiceDialog(parent,'Select file(s) to extract from zip file '+str(filename),
            'Choose file(s)',choices)
        if dlg.ShowModal() == wx.ID_OK:
            zlist = dlg.GetSelections()
        else:
            zlist = []
        dlg.Destroy()
    else:
        # select one from a from list
        choices = [i.filename for i in zinfo]
        dlg = wx.SingleChoiceDialog(parent,
            'Select file to extract from zip file'+str(filename),'Choose file',
            choices,)
        if dlg.ShowModal() == wx.ID_OK:
            zlist = [dlg.GetSelection()]
        else:
            zlist = [-1]
        dlg.Destroy()
        
    outlist = []
    for zindex in zlist:
        if zindex >= 0:
            efil = os.path.join(zloc, os.path.split(zinfo[zindex].filename)[1])
            if os.path.exists(efil) and confirmoverwrite:
                result = wx.ID_NO
                dlg = wx.MessageDialog(parent,
                    'File '+str(efil)+' already exists. OK to overwrite it?',
                    'Confirm overwrite',wx.YES_NO | wx.ICON_QUESTION)
                try:
                    result = dlg.ShowModal()
                finally:
                    dlg.Destroy()
                if result == wx.ID_NO:
                    zindex = -1
        if zindex >= 0:
            # extract the file to the current directory, regardless of it's original path
            #z.extract(zinfo[zindex],zloc)
            eloc,efil = os.path.split(zinfo[zindex].filename)
            outfile = os.path.join(zloc, efil)
            fpin = z.open(zinfo[zindex])
            fpout = open(outfile, "wb")
            shutil.copyfileobj(fpin, fpout)
            fpin.close()
            fpout.close()
            outlist.append(outfile)
    z.close()
    if multipleselect and len(outlist) >= 1:
        return outlist
    elif len(outlist) == 1:
        return outlist[0]
    else:
        return None

def askQuestion(parent,question,title):
    '''Simple code to ask a Y/N question and get answer'''
    ans = True
    try:
        dlg = wx.MessageDialog(parent,question,title,wx.YES_NO | wx.ICON_QUESTION)
        ans = (dlg.ShowModal() == wx.ID_YES)
    finally:
        dlg.Destroy()
    return ans

#===========================================================================
def gitFetch(G2frame):
    wx.BeginBusyCursor()
    pdlg = wx.ProgressDialog('Updating','Performing git update',11,
                    style = wx.PD_ELAPSED_TIME|wx.PD_CAN_ABORT,
                    parent=G2frame)
    try:
        pdlg.CenterOnParent()
        if hasattr(G2frame,'UpdateTask'):     # check if git update has completed
            count = 0
            while G2frame.UpdateTask.poll() is None:
                count += 1
                if count > 10:
                    G2MessageBox(G2frame,
                        'Background git update has not completed, try again later', 
                        title='Warning')
                    return
                time.sleep(1)
                ok,_ = pdlg.Update(count)
                wx.GetApp().Yield()
                if not ok:
                    return
        if GSASIIpath.GetConfigValue('debug'): print('background update complete')
        # try update one more time just to make sure
        GSASIIpath.gitGetUpdate('immediate')
    except Exception as msg:
        raise Exception(msg)
    finally:
        pdlg.Destroy()
        wx.EndBusyCursor()
            
def gitCheckUpdates(G2frame):
    '''Used to update to the latest GSAS-II version, but checks for a variety
    of repository conditions that could make this process more complex. If 
    there are uncommitted local changes, these changes must be cached or 
    deleted first. If there are local changes that have been committed or a new 
    branch has been created, the user (how obstensibly must know use of git)
    will probably need to do this manually. If GSAS-II has previously been 
    regressed (using :func:`gitSelectVersion`), then this is noted as well.

    When all is done, function :func:`GSASIIpath.gitStartUpdate` is called to 
    actually perform the update.
    '''
    try:
        gitFetch(G2frame)  # download latest updates from server
    except:
        G2MessageBox(G2frame,
                    'Unable to access updates: no internet connection?',
                    title='git error')
        return

    status = GSASIIpath.gitTestGSASII()
    if status < 0:
        G2MessageBox(G2frame,
                    'Problem with git access to GSAS-II. Seek help or reinstall',
                    title='git error')
        return
    localChanges = bool(status & 5) # If True need to select between
                                    # --git-stash or --git-reset
                                    # False: neither should be used
    
    if status&8:  # not on local branch
        if localChanges:
            G2MessageBox(G2frame,
                    'You have made a local branch and have uncommited changes. Save, stash or restore then before an update can be done.',
                    title='update not possible')
            return
        else:
            msg = ('You have made a local branch and have switched to that.'+
                   ' Do you want to update anyway? Doing so will return'+
                   ' you to the master branch. This may be better handled'+
                   ' manually.\n\nPress "Yes" to continue with update\n'+
                   'Press "Cancel" to stop the update.')
        dlg = wx.MessageDialog(G2frame, msg, 'Confirm update?',
                wx.OK|wx.CANCEL|wx.CANCEL_DEFAULT|wx.CENTRE|wx.ICON_QUESTION)
        ans = dlg.ShowModal()
        dlg.Destroy()
        if ans == wx.ID_CANCEL: return

    regressmsg = ''
    if status&2:  # detached head
        remoteupdates,localupdates = GSASIIpath.gitCountRegressions()
        if localupdates:
            msg = ("Your copy of GSAS-II has been regressed. You are "+
                    f"{remoteupdates} versions behind the current. "+
                    f"WARNING: you have committed {localupdates} "+
                    "changes locally. Your local commits will be difficult "+
                    "to locate if you update. SUGGESTION: if your changes "+
                    "are meaningful, create a new git branch and return "+
                    "to the master branch.\n\n"+
                    'Press "Yes" to continue, stranding your local changes\n'+
                    'Press "Cancel" to stop the update.')
            dlg = wx.MessageDialog(G2frame, msg, 'Confirm update?',
                    wx.YES|wx.CANCEL|wx.CANCEL_DEFAULT|wx.CENTRE|wx.ICON_QUESTION)
            ans = dlg.ShowModal()
            dlg.Destroy()
            if ans == wx.ID_CANCEL: return
        regressmsg = f"Your copy of GSAS-II has been regressed. You are {remoteupdates} versions behind the current.\n\n"
    else:
        # on head, standard update, is update needed or blocked by local changes
        rc,lc,_ = GSASIIpath.gitCheckForUpdates(False)
        if len(rc) == 0:
            G2MessageBox(G2frame,
                    'Your copy of GSAS-II is up to date. No update is needed.',
                    title='no updates')
            return
        if len(lc) != 0:
            msg = ('You have made local changes and committed them '+
                   'into the master branch. GUI-based updates cannot '+
                   'be made. You should perform a git merge manually')
            G2MessageBox(G2frame,msg,title='Do manual update')
            return

    cmds = ['--git-update']
    if localChanges:
        if gitAskLocalChanges(G2frame,cmds): return
    if gitAskSave(G2frame,regressmsg,cmds): return
    # launch changes and restart
    GSASIIpath.gitStartUpdate(cmds)

def gitAskLocalChanges(G2frame,cmds):
    msg = ('You have locally-made changes to the GSAS-II source '+
            'code files. Do you want to discard these changes?\n\n'+
            'If you select "Yes" the changes will be overwritten '+
            'before the update.\nIf you select "No" the changes will '+
            'be archived (git stash).\nSelect "Cancel" '+
            ' to discontinue the update process.')
    dlg = wx.MessageDialog(G2frame, msg, 'Save/Discard Local Changes?',
            wx.YES_NO|wx.CANCEL|wx.NO_DEFAULT|wx.CENTRE|wx.ICON_QUESTION)
    ans = dlg.ShowModal()
    dlg.Destroy()
    if ans == wx.ID_CANCEL:
        return True
    elif ans == wx.ID_YES:
        cmds += ['--git-reset']
    else:
        cmds += ['--git-stash']
    return False

def gitAskSave(G2frame,regressmsg,cmds):
    # before doing update, need to save project
    msg = ('Before continuing, do you want to save your project? '+
           'Select "Yes" to save, "No" to skip the save, or "Cancel"'+
           ' to discontinue the update process.\n\n'+
           'If "Yes", GSAS-II will reopen the project after the update. '+
           'The update will now begin unless Cancel is pressed.')
    dlg = wx.MessageDialog(G2frame, regressmsg+msg,
                'Save Project and Start Update?',
                wx.YES_NO|wx.CANCEL|wx.YES_DEFAULT|wx.CENTRE|wx.ICON_QUESTION)
    ans = dlg.ShowModal()
    dlg.Destroy()
    if ans == wx.ID_CANCEL:
        return True
    elif ans == wx.ID_YES:
        ans = G2frame.OnFileSave(None)
        if ans: 
            cmds += [G2frame.GSASprojectfile]
        return False
    
def gitSelectVersion(G2frame):
    '''Used to regress to a previous GSAS-II version, checking first 
    for a variety of repository conditions that could make this process 
    more complex. If there are uncommitted local changes, these changes 
    must be cached or deleted before a different version can be installed. 
    If there are local changes that have been committed or a new 
    branch has been created, the user (how obstensibly must know use of git)
    will probably need to do this manually. If GSAS-II has previously been 
    regressed (using :func:`gitSelectVersion`), then this is noted as well.

    When all is done, function :func:`GSASIIpath.gitStartUpdate` is called to 
    actually perform the update.
    '''
    # get updates from server
    gitFetch(G2frame)  # download latest updates from server

    status = GSASIIpath.gitTestGSASII()
    if status < 0:
        G2MessageBox(G2frame,
                    'Problem with git access to GSAS-II. Seek help or reinstall',
                    title='git error')
        return
    localChanges = bool(status & 5) # If True need to select between
                                    # --git-stash or --git-reset
                                    # False: neither should be used
    if status&8:  # not on local branch
        if localChanges:
            msg =  ('You have switched to a local branch and have '+
                    'uncommited changes. Save, stash or restore and switch '+
                    'back to the master branch. Regression is not possible '+
                    'from the GUI when on any other branch.')
        else:
            msg =  ('You have made and switched to a local branch. '+
                    'You must manually switch back to the master branch.'+
                    ' Regression is not possible '+
                    'from the GUI when on any other branch.')
        G2MessageBox(G2frame,msg,title='update not possible')
        return

    regressmsg = ''
    if status&2:  # detached head
        remoteupdates,localupdates = GSASIIpath.gitCountRegressions()
        if localupdates:
            msg = ("Your copy of GSAS-II has been regressed. You are "+
                    f"{remoteupdates} versions behind the current. "+
                    f"WARNING: you have committed {localupdates} "+
                    "changes locally. Your local commits will be difficult "+
                    "to locate if you update. SUGGESTION: if your changes "+
                    "are meaningful, create a new git branch and return "+
                    "to the master branch.\n\n"+
                    'Press "Yes" to continue, stranding your local changes\n'+
                    'Press "Cancel" to stop the update.')
            dlg = wx.MessageDialog(G2frame, msg, 'Confirm update?',
                    wx.YES|wx.CANCEL|wx.CANCEL_DEFAULT|wx.CENTRE|wx.ICON_QUESTION)
            ans = dlg.ShowModal()
            dlg.Destroy()
            if ans == wx.ID_CANCEL: return
        regressmsg = f"Your copy of GSAS-II has been regressed. You are {remoteupdates} versions behind the current.\n\n"
    else:
        # on head, standard update, is update blocked by local changes
        rc,lc,_ = GSASIIpath.gitCheckForUpdates(False)
        if len(lc) != 0:
            msg = ('You have made local changes and committed them '+
                   'into the master branch. GUI-based regression cannot '+
                   'be done. You should perform a "git checkout" to the '+
                   'desired version manually')
            G2MessageBox(G2frame,msg,title='Do manual update')
            return
        
    # browse and select a version here
    dlg = gitVersionSelector()
    ans = dlg.ShowModal()
    if ans == wx.ID_CANCEL: return
    githash = dlg.getVersion()
    if githash is None:
        print('Nothing to be done')
        return
    if githash == 0:
        print('select newest GSAS-II version')
        cmds = ['--git-update']
    else:
        cmds = [f'--git-regress={githash}']

    if localChanges:
        if gitAskLocalChanges(G2frame,cmds): return
    if gitAskSave(G2frame,regressmsg,cmds): return

    # launch changes and restart
    GSASIIpath.gitStartUpdate(cmds)

def gitSelectBranch(event):
    '''Pull in latest GSAS-II branches on origin server; Allow user to 
    select a branch; checkout that branch and restart GSAS-II. 
    Expected to be used by developers and by expert users only.
    '''
    G2frame = wx.App.GetMainTopWindow()
    gitInst = GSASIIpath.HowIsG2Installed()
    if not gitInst.startswith('github-rev'):
        G2MessageBox(G2frame,
            'Unable to switch branches unless GSAS-II has been installed from GitHub; installed as: '+gitInst,
            'Not a git install')
        return
    if not os.path.exists(GSASIIpath.path2GSAS2): 
        print(f'Warning: Directory {GSASIIpath.path2GSAS2} not found')
        return
    if os.path.exists(os.path.join(GSASIIpath.path2GSAS2,'..','.git')):
        path2repo = os.path.join(path2GSAS2,'..')  # expected location
    elif os.path.exists(os.path.join(GSASIIpath.path2GSAS2,'.git')):
        path2repo = GSASIIpath.path2GSAS2
    else:
        print(f'Warning: Repository {path2GSAS2} not found')
        return
    try:
        g2repo = GSASIIpath.openGitRepo(path2repo)
    except Exception as msg:
        print(f'Warning: Failed to open repository. Error: {msg}')
        return
    if g2repo.is_dirty() or g2repo.index.diff("HEAD"): # changed or staged files
        G2MessageBox(G2frame,
            'You have local changes. They must be reset, committed or stashed before switching branches',
            'Local changes')
        return
    if g2repo.head.is_detached:
        G2MessageBox(G2frame,
            'You have a old previous version loaded; you must be on a branch head to switching branches',
            'Detached head')
        return

    # make sure that branches are accessible & get updates
    print('getting updates...',end='')
    g2repo.git.remote('set-branches','origin','*')
    print('..',end='')
    g2repo.git.fetch()
    print('.done')
    branchlist = [i.strip() for i in g2repo.git.branch('-r').split('\n') if '->' not in i]
    choices = [i for i in  [os.path.split(i)[1] for i in branchlist] if i != g2repo.active_branch.name]
    if len(choices) == 0:
        G2MessageBox(G2frame,
            'No branches were found to select. Unexpected!',
            'No branches')
        return
    if len(choices) == 1:
        b = choices[0]
    else:
        dlg = G2SingleChoiceDialog(G2frame,'Select branch to use','Select Branch',
                                 choices)
        dlg.CenterOnParent()
        try:
            if dlg.ShowModal() == wx.ID_OK:
                b = choices[dlg.GetSelection()]
            else:
                return
        finally:
            dlg.Destroy()
    msg = f'''Confirm switching from git branch {g2repo.active_branch.name!r} to {b!r}.

If confirmed here, GSAS-II will restart. 

Do you want to save your project before restarting?
Select "Yes" to save, "No" to skip the save, or "Cancel"
to discontinue the restart process.

If "Yes", GSAS-II will reopen the project after the update.

The switch will be made unless Cancel is pressed.'''
    dlg = wx.MessageDialog(G2frame, msg, 'Confirm branch switch?',
                wx.YES_NO|wx.CANCEL|wx.YES_DEFAULT|wx.CENTRE|wx.ICON_QUESTION)
    ans = dlg.ShowModal()
    dlg.Destroy()
    if ans == wx.ID_CANCEL:
        return
    elif ans == wx.ID_YES:
        ans = G2frame.OnFileSave(None)
        if not ans: return
        project = os.path.abspath(G2frame.GSASprojectfile)
        print(f"Restarting GSAS-II with project file {project!r}")
    else:
        print("Restarting GSAS-II without a project file ")
        project = None
    # I hope that it is possible to do a checkout on Windows
    # (source files are not locked). If this is not the case
    # then another approach will be needed, where a .bat file is used
    # or GSASIIpath is used, as is the case for updates
    a = g2repo.git.checkout(b)
    if 'Your branch is behind' in a:
        print('updating local copy of branch')
        print(g2repo.git.pull())
    G2fil.openInNewTerm(project)
    print ('exiting GSAS-II')
    sys.exit()
    
#===========================================================================
def svnCheckUpdates(G2frame):
    '''Check if the GSAS-II repository has an update for the current 
    source files and perform that update if requested.
    '''
    wx.BeginBusyCursor()
    local = GSASIIpath.svnGetRev()
    if local is None: 
        wx.EndBusyCursor()
        dlg = wx.MessageDialog(G2frame,
                               'Unable to run subversion on the GSAS-II current directory. Is GSAS-II installed correctly?',
                               'Subversion error',
                               wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
        return
    print ('Installed GSAS-II version: '+local)
    repos = GSASIIpath.svnGetRev(local=False)
    wx.EndBusyCursor()
    # has the current branch disappeared? If so, switch to the trunk -- not fully tested
    if (repos is None and "not found" in GSASIIpath.svnLastError.lower()
        and "path" in GSASIIpath.svnLastError.lower()):
        print('Repository is gone, will switch to trunk')
        GSASIIpath.svnSwitch2branch()
        return
    elif repos is None: 
        dlg = wx.MessageDialog(G2frame,
                               'Unable to access the GSAS-II server. Is this computer on the internet?',
                               'Server unavailable',
                               wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
        return
    errmsg,warnmsg = G2gd.TestOldVersions()
    if (errmsg or warnmsg):
        msg = 'Based on the age of Python or an installed Python package (see below)'
        msg += ' you are recommended to'
        if GSASIIpath.condaTest():
            msg += ' either use conda to update (see https://bit.ly/G2pkgs for version recommendations) and then update GSAS-II. Or'
        msg += ' reinstall GSAS-II (see https://bit.ly/G2install), which will update both.\n\n'
        if errmsg:
            opt = wx.YES_NO|wx.ICON_QUESTION|wx.CANCEL|wx.NO_DEFAULT
            msg += 'Error(s):\n\t'+errmsg
        else:
            opt = wx.YES_NO|wx.ICON_QUESTION|wx.CANCEL|wx.YES_DEFAULT
        if warnmsg:
            if errmsg: msg += '\n\nWarning(s):\n\t'
            msg += warnmsg

        msg += '\n\nContinue to update GSAS-II?'
        dlg = wx.MessageDialog(G2frame, msg,'Confirm update',opt)
        result = wx.ID_NO
        try:
            result = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if result != wx.ID_YES: return
    print ('GSAS-II version on server: '+repos)
    if local == repos:
        up,mod,lock = GSASIIpath.svnGetFileStatus()
    else:
        up = 0
        mods = GSASIIpath.svnFindLocalChanges()
        mod = len(mods)
    if local == repos and up:
        dlg = wx.MessageDialog(G2frame,
                               'You have the current version '+
                               ' of GSAS-II installed ('+repos+
                               '). However, '+str(up)+
                               ' file(s) still need to be updated.'+
                               ' Most likely a previous update failed to finish.'+
                               '\n\nPress OK to retry the update.'+
                               '\n\nIf this problem continues contact toby@anl.gov',
                               'Failed Update Likely',
                               wx.OK|wx.CANCEL)
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        else:
            dlg.Destroy()
        if lock: GSASIIpath.svnCleanup()
    elif local == repos:
        dlg = wx.MessageDialog(G2frame,
                               'GSAS-II is up-to-date. Version '+local+' is already loaded.',
                               'GSAS-II Up-to-date',
                               wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
        return
    if mod:
        dlg = wx.MessageDialog(G2frame,
                               'You have version '+local+
                               ' of GSAS-II installed, but the current version is '+repos+
                               '. However, '+str(mod)+
                               ' file(s) on your local computer have been modified.'
                               ' Updating will attempt to merge your local changes with '
                               'the latest GSAS-II version, but if '
                               'conflicts arise, local changes will be '
                               'discarded. It is also possible that the '
                               'merge may prevent GSAS-II from running. '
                               '\n\nPress OK to start an update if this is acceptable:',
                               'Local GSAS-II Mods',
                               wx.OK|wx.CANCEL)
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        else:
            dlg.Destroy()
    else:
        dlg = wx.MessageDialog(G2frame,
                               'You have version '+local+
                               ' of GSAS-II installed, but the current version is '+repos+
                               '. Press OK to start an update:',
                               'GSAS-II Updates',
                               wx.OK|wx.CANCEL)
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        dlg.Destroy()
    print ('start updates')
    if G2frame.GPXtree.GetCount() > 1:
        dlg = wx.MessageDialog(G2frame,
                           'Your project will now be saved, GSAS-II will exit and an update '+
                           'will be performed and GSAS-II will restart. Press Cancel '+
                           'in next dialog to avoid saving the project.',
                           'Starting update',
                           wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
        G2frame.OnFileSave(None)
        GPX = G2frame.GSASprojectfile
        GSASIIpath.svnUpdateProcess(projectfile=GPX)
    else:
        GSASIIpath.svnUpdateProcess()
    return

def svnSelectVersion(G2frame):
    '''Allow the user to select a specific version of GSAS-II from the 
    APS svn server
    '''    
    local = GSASIIpath.svnGetRev()
    if local is None:
        dlg = wx.MessageDialog(G2frame,
                               'Unable to run subversion on the GSAS-II current directory. Is GSAS-II installed correctly?',
                               'Subversion error',
                               wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
        return
    mods = GSASIIpath.svnFindLocalChanges()
    if mods:
        dlg = wx.MessageDialog(G2frame,
                               'You have version '+local+
                               ' of GSAS-II installed'
                               '. However, '+str(len(mods))+
                               ' file(s) on your local computer have been modified.'
                               ' Downdating will attempt to merge your local changes with '
                               'the selected GSAS-II version. '
                               'Downdating is not encouraged because '
                               'if merging is not possible, your local changes will be '
                               'discarded. It is also possible that the '
                               'merge may prevent GSAS-II from running. '
                               'Press OK to continue anyway.',
                               'Local GSAS-II Mods',
                               wx.OK|wx.CANCEL)
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        dlg.Destroy()
    if GSASIIpath.svnGetRev(local=False) is None:
        dlg = wx.MessageDialog(G2frame,
                               'Error obtaining current GSAS-II version. Is internet access working correctly?',
                               'Subversion error',
                               wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
        return
    dlg = downdate(parent=G2frame)
    if dlg.ShowModal() == wx.ID_OK:
        ver = dlg.getVersion()
    else:
        dlg.Destroy()
        return
    dlg.Destroy()
    print('start regress to '+str(ver))
    G2frame.OnFileSave(None)
    GPX = G2frame.GSASprojectfile
    GSASIIpath.svnUpdateProcess(projectfile=GPX,version=str(ver))
    return

# Importer GUI stuff
def ImportMsg(parent,msgs):
    '''Show a message with the warnings from importers that 
    could not be installed (due to uninstalled Python packages). Then 
    offer the chance to install GSAS-II packages using :func:`SelectPkgInstall`
    '''
    text = ('Message(s) from load of importers\n\n  '+
                '\n\n'.join(msgs)+
                '\n\nNote: These errors only need to be addressed if you want to use the importers listed above')
    ShowScrolledInfo(parent,text,
                    header='Importer load problems',
                    width=650,
                    buttonlist=[('Install packages',SelectPkgInstall), wx.ID_CLOSE]
                         )

def SelectPkgInstall(event):
    '''Offer the user a chance to install Python packages needed by one or 
    more importers. There might be times where something like this will be 
    useful for other GSAS-II actions.
    '''
    dlg = event.GetEventObject().GetParent()
    dlg.EndModal(wx.ID_OK)
    G2frame = wx.App.GetMainTopWindow()
    choices = []
    for key in G2fil.condaRequestList:
        for item in G2fil.condaRequestList[key]:
            choices.append((item,key))
    msg = 'Select packages to install'
    if GSASIIpath.condaTest():
        msg += ' using conda'
    else:
        msg += ' using pip'
    sel = MultiColMultiSelDlg(G2frame, 'Install packages?', msg,
                             [('package',120,0),('needed by',300,0)], choices)
    if sel is None: return
    if not any(sel): return
    pkgs = [choices[i][0] for i,f in enumerate(sel) if f]
    if GSASIIpath.condaTest():
        if not GSASIIpath.condaTest(True):
            GSASIIpath.addCondaPkg()
        err = GSASIIpath.condaInstall(pkgs)
        if err:
            print(f'Error from conda: {err}')
            return
    else:
        err = GSASIIpath.pipInstall(pkgs)
        if err:
            print(f'Error from pip: {err}')
            return
    msg = '''You must restart GSAS-II to access the importer(s) 
requiring the installed package(s). 

Select "Yes" to save, "No" to skip the save, or "Cancel"
to discontinue the restart process and continue GSAS-II 
without the importer(s). 

If "Yes", GSAS-II will reopen the project after the update.
'''
    dlg = wx.MessageDialog(G2frame, msg, 'Save and restart?',
                wx.YES_NO|wx.CANCEL|wx.YES_DEFAULT|wx.CENTRE|wx.ICON_QUESTION)
    ans = dlg.ShowModal()
    dlg.Destroy()
    if ans == wx.ID_CANCEL:
        return
    elif ans == wx.ID_YES:
        ans = G2frame.OnFileSave(None)
        if not ans: return
        project = os.path.abspath(G2frame.GSASprojectfile)
        print(f"Restarting GSAS-II with project file {project!r}")
    else:
        print("Restarting GSAS-II without a project file ")
        project = None
    G2fil.openInNewTerm(project)
    print ('exiting GSAS-II')
    sys.exit()

if __name__ == '__main__':
    app = wx.App()
    GSASIIpath.InvokeDebugOpts()
    frm = wx.Frame(None) # create a frame
    ms = wx.BoxSizer(wx.VERTICAL)
    #siz = G2SliderWidget(pnl,valArr,'k','test slider w/entry',.2,1.2,100)
    #ms.Add(siz)
    #siz = G2SliderWidget(pnl,valArr,'k','test slider w/entry',20,50,.1)
    #ms.Add(siz)
    text = 'this is a long string that will be scrolled'
    ms.Add(ScrolledStaticText(frm,label=text))
    txt = ScrolledStaticText(frm,label=text, lbllen=20)
    smallfont = wx.SystemSettings.GetFont(wx.SYS_SYSTEM_FONT)
    smallfont.SetPointSize(10)
    txt.SetFont(smallfont)
    ms.Add(txt)
    ms.Add(ScrolledStaticText(frm,label=text,dots=False,delay=250, lbllen=20))
    frm.SetSizer(ms)
    frm.Show(True)

    G2frame = frm


    testAtoms = ['']
    
    nm = [' ','0','1','-1','2','-2','3','-3','4','5','6','7','8','9']
    dm = ['1','2','3','4','5','6']
    kfmt = ['choice','/','choice',',    ','choice','/','choice',',    ','choice','/','choice',' ']
    def strTest(text):
            if '.' in text: # no decimals
                return False
            elif text.strip() in  [' ','0','1','-1','3/2']: # specials
                return True
            elif '/' in text: #process fraction 
                nums = text.split('/')
                return (0 < int(nums[1]) < 10) and (0 < abs(int(nums[0])) < int(nums[1]))
            return False

    msg = 'test of MultiDataDialog'
    kvec = [['0','0','0'],[' ',' ',' '],[' ',' ',' ',' ']]
    dlg = MultiDataDialog(G2frame,title='k-SUBGROUPSMAG options',
            prompts=[' k-vector 1',' k-vector 2',' k-vector 3',
                     ' Use whole star',' Filter by','preserve axes',
                     'test for mag. atoms','all have moment','max unique'],
            values=kvec+[False,'',True,'',False,100],
            limits=[['0','0','0'],['0','0','0'],['0','0','0'],
                    [True,False],['',' Landau transition',' Only maximal subgroups',],
                [True,False],testAtoms,[True,False],[1,100]],
            testfxns = [[strTest,strTest,strTest],
                        [strTest,strTest,strTest],
                        [strTest,strTest,strTest],
                        None,None,None,None,None,None],
            formats=[['testfxn','testfxn','testfxn'],
                     ['testfxn','testfxn','testfxn'],
                     ['testfxn','testfxn','testfxn'],
                     'bool','choice','bool','choice','bool','%d',],
            header=msg)
    if dlg.ShowModal() == wx.ID_OK: print(dlg.GetValues())
    

    # if True:
    #   title='title here'
    #   header = 'this is where a header goes. this is where a header to explain what to do goes this is where a header to explain what to do goes'
    #   choices = [('xmltodict', 'Bruker .brml Importer'),
    #              ('zarr', 'MIDAS Zarr importer'),
    #              ('h5py', 'HDF5 image importer'),
    #              ('hdf5', 'HDF5 image importer'),
    #              ('test',),('val','used','ignored'),[]]
    #   colInfo = [('package',100, False),
    #                ('needed by',200, False)]

    #   parent = wx.App.GetMainTopWindow()
    #   print(MultiColMultiSelDlg(parent, title, header, colInfo, choices))

    sys.exit()
    app.MainLoop()
    
#    choices = [wx.ID_YES,wx.ID_NO]
#    warnmsg = '\nsome info\non a few lines\nand one more'
#    ans = ShowScrolledInfo(header='Constraint Warning',
#                    txt='Warning noted after last constraint edit:\n'+warnmsg+
#                    '\n\nKeep this change?',
#                    buttonlist=choices,parent=frm,height=250)
#    print(ans, choices)

    import sys; sys.exit()
    
    #======================================================================
    # test Grid with GridFractionEditor
    #======================================================================
    # tbl = [[1.,2.,3.],[1.1,2.1,3.1]]
    # colTypes = 3*[wg.GRID_VALUE_FLOAT+':10,5',]
    # Gtbl = Table(tbl,types=colTypes,rowLabels=['a','b'],colLabels=['1','2','3'])
    # Grid = GSGrid(frm)
    # Grid.SetTable(Gtbl,True)
    # for i in (0,1,2):
    #     attr = wx.grid.GridCellAttr()
    #     attr.IncRef()
    #     attr.SetEditor(GridFractionEditor(Grid))
    #     Grid.SetColAttr(i, attr)
    # frm.SetSize((400,200))
    # app.MainLoop()
    # sys.exit()
    #======================================================================
    # test Tutorial access
    #======================================================================
    # dlg = OpenTutorial(frm)
    # if dlg.ShowModal() == wx.ID_OK:
    #     print("OK")
    # else:
    #     print("Cancel")
    # dlg.Destroy()
    # sys.exit()
    #======================================================================
    # test ScrolledMultiEditor
    #======================================================================
    # Data1 = {
    #      'Order':1,
    #      'omega':'string',
    #      'chi':2.0,
    #      'phi':'',
    #      }
    # elemlst = sorted(Data1.keys())
    # prelbl = sorted(Data1.keys())
    # dictlst = len(elemlst)*[Data1,]
    #Data2 = [True,False,False,True]
    #Checkdictlst = len(Data2)*[Data2,]
    #Checkelemlst = range(len(Checkdictlst))
    # print 'before',Data1,'\n',Data2
    # dlg = ScrolledMultiEditor(
    #     frm,dictlst,elemlst,prelbl,
    #     checkdictlst=Checkdictlst,checkelemlst=Checkelemlst,
    #     checklabel="Refine?",
    #     header="test")
    # if dlg.ShowModal() == wx.ID_OK:
    #     print "OK"
    # else:
    #     print "Cancel"
    # print 'after',Data1,'\n',Data2
    # dlg.Destroy()
    # Data3 = {
    #      'Order':1.0,
    #      'omega':1.1,
    #      'chi':2.0,
    #      'phi':2.3,
    #      'Order1':1.0,
    #      'omega1':1.1,
    #      'chi1':2.0,
    #      'phi1':2.3,
    #      'Order2':1.0,
    #      'omega2':1.1,
    #      'chi2':2.0,
    #      'phi2':2.3,
    #      }
    # elemlst = sorted(Data3.keys())
    # dictlst = len(elemlst)*[Data3,]
    # prelbl = elemlst[:]
    # prelbl[0]="this is a much longer label to stretch things out"
    # Data2 = len(elemlst)*[False,]
    # Data2[1] = Data2[3] = True
    # Checkdictlst = len(elemlst)*[Data2,]
    # Checkelemlst = range(len(Checkdictlst))
    #print 'before',Data3,'\n',Data2
    #print dictlst,"\n",elemlst
    #print Checkdictlst,"\n",Checkelemlst
    # dlg = ScrolledMultiEditor(
    #     frm,dictlst,elemlst,prelbl,
    #     checkdictlst=Checkdictlst,checkelemlst=Checkelemlst,
    #     checklabel="Refine?",
    #     header="test",CopyButton=True)
    # if dlg.ShowModal() == wx.ID_OK:
    #     print "OK"
    # else:
    #     print "Cancel"
    #print 'after',Data3,'\n',Data2

    # Data2 = list(range(100))
    # elemlst += range(2,6)
    # postlbl += range(2,6)
    # dictlst += len(range(2,6))*[Data2,]

    # prelbl = range(len(elemlst))
    # postlbl[1] = "a very long label for the 2nd item to force a horiz. scrollbar"
    # header="""This is a longer\nmultiline and perhaps silly header"""
    # dlg = ScrolledMultiEditor(frm,dictlst,elemlst,prelbl,postlbl,
    #                           header=header,CopyButton=True)
    # print Data1
    # if dlg.ShowModal() == wx.ID_OK:
    #     for d,k in zip(dictlst,elemlst):
    #         print k,d[k]
    # dlg.Destroy()
    # if CallScrolledMultiEditor(frm,dictlst,elemlst,prelbl,postlbl,
    #                            header=header):
    #     for d,k in zip(dictlst,elemlst):
    #         print k,d[k]

    #======================================================================
    # test G2MultiChoiceDialog
    #======================================================================
    # choices = []
    # for i in range(21):
    #     choices.append("option_"+str(i))
    # od = {
    #     'label_1':'This is a bool','value_1':True,
    #     'label_2':'This is a int','value_2':-1,
    #     'label_3':'This is a float','value_3':1.0,
    #     'label_4':'This is a string','value_4':'test',}
    # dlg = G2MultiChoiceDialog(frm, 'Sequential refinement',
    #                           'Select dataset to include',
    #                           choices,extraOpts=od)
    # sel = range(2,11,2)
    # dlg.SetSelections(sel)
    # dlg.SetSelections((1,5))
    # if dlg.ShowModal() == wx.ID_OK:
    #     for sel in dlg.GetSelections():
    #         print (sel,choices[sel])
    # print (od)
    # od = {}
    # dlg = G2MultiChoiceDialog(frm, 'Sequential refinement',
    #                           'Select dataset to include',
    #                           choices,extraOpts=od)
    # sel = range(2,11,2)
    # dlg.SetSelections(sel)
    # dlg.SetSelections((1,5))
    # if dlg.ShowModal() == wx.ID_OK: pass
    #======================================================================
    # test wx.MultiChoiceDialog
    #======================================================================
    # dlg = wx.MultiChoiceDialog(frm, 'Sequential refinement',
    #                           'Select dataset to include',
    #                           choices)
    # sel = range(2,11,2)
    # dlg.SetSelections(sel)
    # dlg.SetSelections((1,5))
    # if dlg.ShowModal() == wx.ID_OK:
    #     for sel in dlg.GetSelections():
    #         print sel,choices[sel]

    # pnl = wx.Panel(frm)
    # siz = wx.BoxSizer(wx.VERTICAL)
    # td = {'Goni':200.,'a':1.,'int':1,'calc':1./3.,'string':'s'}
    # for key in sorted(td):
    #     txt = ValidatedTxtCtrl(pnl,td,key,typeHint=float)
    #     siz.Add(txt)
    # pnl.SetSizer(siz)
    # siz.Fit(frm)
    # app.MainLoop()
    # print(td)
    choicelist=[ ('a','b','c'), ('test1','test2'),('no choice',)]
    headinglist = [ 'select a, b or c', 'select 1 of 2', 'No option here']
    dlg = MultipleChoicesDialog(choicelist,headinglist,parent=frm)
    if dlg.ShowModal() == wx.ID_OK:
        print(dlg.chosen)
    print(MultipleChoicesSelector(choicelist,headinglist,frm))
    pnl = wx.Panel(frm)
    valArr = {'k':1.0}
    ms = wx.BoxSizer(wx.VERTICAL)
    #siz = G2SliderWidget(pnl,valArr,'k','test slider w/entry',.2,1.2,100)
    #ms.Add(siz)
    siz = G2SliderWidget(pnl,valArr,'k','test slider w/entry',2,5,1)
    ms.Add(siz)
    #siz = G2SliderWidget(pnl,valArr,'k','test slider w/entry',20,50,.1)
    #ms.Add(siz)
    pnl.SetSizer(ms)
    ms.Fit(frm)
    app.MainLoop()
    print(valArr)
