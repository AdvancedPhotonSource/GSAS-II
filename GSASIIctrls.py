# -*- coding: utf-8 -*-
#GSASIIctrls - Custom GSAS-II GUI controls
########### SVN repository information ###################
# $Date: $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
*GSASIIctrls: Custom GUI controls*
-------------------------------------------

A library of GUI controls for reuse throughout GSAS-II

(at present many are still in GSASIIgrid, but with time will be moved here)

'''
import wx
# import wx.grid as wg
# import wx.wizard as wz
# import wx.aui
import wx.lib.scrolledpanel as wxscroll
import time
import copy
# import cPickle
import sys
import os
# import numpy as np
# import numpy.ma as ma
# import scipy.optimize as so
import wx.html        # could postpone this for quicker startup
import webbrowser     # could postpone this for quicker startup

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1614 $")
# import GSASIImath as G2mth
# import GSASIIIO as G2IO
# import GSASIIstrIO as G2stIO
# import GSASIIlattice as G2lat
# import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
# import GSASIIimgGUI as G2imG
# import GSASIIphsGUI as G2phG
# import GSASIIspc as G2spc
# import GSASIImapvars as G2mv
# import GSASIIconstrGUI as G2cnstG
# import GSASIIrestrGUI as G2restG
import GSASIIpy3 as G2py3
# import GSASIIobj as G2obj
# import GSASIIexprGUI as G2exG
import GSASIIlog as log

# Define a short names for convenience
WHITE = (255,255,255)
DULL_YELLOW = (230,230,190)
VERY_LIGHT_GREY = wx.Colour(235,235,235)
WACV = wx.ALIGN_CENTER_VERTICAL

################################################################################
#### Tree Control
################################################################################
class G2TreeCtrl(wx.TreeCtrl):
    '''Create a wrapper around the standard TreeCtrl so we can "wrap"
    various events.
    
    This logs when a tree item is selected (in :meth:`onSelectionChanged`)

    This also wraps lists and dicts pulled out of the tree to track where
    they were retrieved from.
    '''
    def __init__(self,parent=None,*args,**kwargs):
        super(self.__class__,self).__init__(parent=parent,*args,**kwargs)
        self.G2frame = parent.GetParent()
        self.root = self.AddRoot('Loaded Data: ')
        self.SelectionChanged = None
        log.LogInfo['Tree'] = self

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

    def onSelectionChanged(self,event):
        '''Log each press on a tree item here. 
        '''
        if self.SelectionChanged:
            textlist = self._getTreeItemsList(event.GetItem())
            if log.LogInfo['Logging'] and event.GetItem() != self.root:
                textlist[0] = self.GetRelativeHistNum(textlist[0])
                if textlist[0] == "Phases" and len(textlist) > 1:
                    textlist[1] = self.GetRelativePhaseNum(textlist[1])
                log.MakeTreeLog(textlist)
            self.SelectionChanged(event)

    def Bind(self,eventtype,handler,*args,**kwargs):
        '''Override the Bind() function so that page change events can be trapped
        '''
        if eventtype == wx.EVT_TREE_SEL_CHANGED:
            self.SelectionChanged = handler
            wx.TreeCtrl.Bind(self,eventtype,self.onSelectionChanged)
            return
        wx.TreeCtrl.Bind(self,eventtype,handler,*args,**kwargs)

    # commented out, disables Logging
    # def GetItemPyData(self,*args,**kwargs):
    #    '''Override the standard method to wrap the contents
    #    so that the source can be logged when changed
    #    '''
    #    data = super(self.__class__,self).GetItemPyData(*args,**kwargs)
    #    textlist = self._getTreeItemsList(args[0])
    #    if type(data) is dict:
    #        return log.dictLogged(data,textlist)
    #    if type(data) is list:
    #        return log.listLogged(data,textlist)
    #    if type(data) is tuple: #N.B. tuples get converted to lists
    #        return log.listLogged(list(data),textlist)
    #    return data

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

################################################################################
#### TextCtrl that stores input as entered with optional validation
################################################################################
class ValidatedTxtCtrl(wx.TextCtrl):
    '''Create a TextCtrl widget that uses a validator to prevent the
    entry of inappropriate characters and changes color to highlight
    when invalid input is supplied. As valid values are typed,
    they are placed into the dict or list where the initial value
    came from. The type of the initial value must be int,
    float or str or None (see :obj:`key` and :obj:`typeHint`);
    this type (or the one in :obj:`typeHint`) is preserved.

    Float values can be entered in the TextCtrl as numbers or also
    as algebraic expressions using operators + - / \* () and \*\*,
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
      
    :param list nDig: number of digits & places ([nDig,nPlc]) after decimal to use
      for display of float. Alternately, None can be specified which causes
      numbers to be displayed with approximately 5 significant figures
      (Default=None).

    :param bool notBlank: if True (default) blank values are invalid
      for str inputs.
      
    :param number min: minimum allowed valid value. If None (default) the
      lower limit is unbounded.

    :param number max: maximum allowed valid value. If None (default) the
      upper limit is unbounded

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
         * tc:      (*wx.TextCtrl*)  the TextCtrl name

      The number of keyword arguments may be increased in the future should needs arise,
      so it is best to code these functions with a \*\*kwargs argument so they will
      continue to run without errors

      The default for OnLeave is None, which indicates no function should
      be called.

    :param type typeHint: the value of typeHint is overrides the initial value
      for the dict/list element ``loc[key]``, if set to 
      int or float, which specifies the type for input to the TextCtrl.
      Defaults as None, which is ignored.

    :param bool CIFinput: for str input, indicates that only printable
      ASCII characters may be entered into the TextCtrl. Forces output
      to be ASCII rather than Unicode. For float and int input, allows
      use of a single '?' or '.' character as valid input.

    :param dict OnLeaveArgs: a dict with keyword args that are passed to
      the :attr:`OnLeave` function. Defaults to ``{}``

    :param (other): other optional keyword parameters for the
      wx.TextCtrl widget such as size or style may be specified.

    '''
    def __init__(self,parent,loc,key,nDig=None,notBlank=True,min=None,max=None,
                 OKcontrol=None,OnLeave=None,typeHint=None,
                 CIFinput=False, OnLeaveArgs={}, **kw):
        # save passed values needed outside __init__
        self.result = loc
        self.key = key
        self.nDig = nDig
        self.OKcontrol=OKcontrol
        self.OnLeave = OnLeave
        self.OnLeaveArgs = OnLeaveArgs
        self.CIFinput = CIFinput
        self.type = str
        # initialization
        self.invalid = False   # indicates if the control has invalid contents
        self.evaluated = False # set to True when the validator recognizes an expression
        val = loc[key]
        if isinstance(val,int) or typeHint is int:
            self.type = int
            wx.TextCtrl.__init__(
                self,parent,wx.ID_ANY,
                validator=NumberValidator(int,result=loc,key=key,
                                          min=min,max=max,
                                          OKcontrol=OKcontrol,
                                          CIFinput=CIFinput),
                **kw)
            if val is not None:
                self._setValue(val)
            else: # no default is invalid for a number
                self.invalid = True
                self._IndicateValidity()

        elif isinstance(val,float) or typeHint is float:
            self.type = float
            wx.TextCtrl.__init__(
                self,parent,wx.ID_ANY,
                validator=NumberValidator(float,result=loc,key=key,
                                          min=min,max=max,
                                          OKcontrol=OKcontrol,
                                          CIFinput=CIFinput),
                **kw)
            if val is not None:
                self._setValue(val)
            else:
                self.invalid = True
                self._IndicateValidity()

        elif isinstance(val,str) or isinstance(val,unicode):
            if self.CIFinput:
                wx.TextCtrl.__init__(
                    self,parent,wx.ID_ANY,val,
                    validator=ASCIIValidator(result=loc,key=key),
                    **kw)
            else:
                wx.TextCtrl.__init__(self,parent,wx.ID_ANY,val,**kw)
            if notBlank:
                self.Bind(wx.EVT_CHAR,self._onStringKey)
                self.ShowStringValidity() # test if valid input
            else:
                self.invalid = False
                self.Bind(wx.EVT_CHAR,self._GetStringValue)
        elif val is None:
            raise Exception,("ValidatedTxtCtrl error: value of "+str(key)+
                             " element is None and typeHint not defined as int or float")
        else:
            raise Exception,("ValidatedTxtCtrl error: Unknown element ("+str(key)+
                             ") type: "+str(type(val)))
        # When the mouse is moved away or the widget loses focus,
        # display the last saved value, if an expression
        #self.Bind(wx.EVT_LEAVE_WINDOW, self._onLeaveWindow)
        self.Bind(wx.EVT_TEXT_ENTER, self._onLoseFocus)
        self.Bind(wx.EVT_KILL_FOCUS, self._onLoseFocus)
        # patch for wx 2.9 on Mac
        i,j= wx.__version__.split('.')[0:2]
        if int(i)+int(j)/10. > 2.8 and 'wxOSX' in wx.PlatformInfo:
            self.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)

    def SetValue(self,val):
        if self.result is not None: # note that this bypasses formatting
            self.result[self.key] = val
            log.LogVarChange(self.result,self.key)
        self._setValue(val)

    def _setValue(self,val):
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
            wx.TextCtrl.SetValue(self,str(val))
        elif self.type is float:
            try:
                val = float(val) # convert strings, if needed
            except:
                if self.CIFinput and (val == '?' or val == '.'):
                    pass
                else:
                    self.invalid = True
            if self.nDig:
                wx.TextCtrl.SetValue(self,str(G2py3.FormatValue(val,self.nDig)))
            else:
                wx.TextCtrl.SetValue(self,str(G2py3.FormatSigFigs(val)).rstrip('0'))
        else:
            wx.TextCtrl.SetValue(self,str(val))
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
        if key == wx.WXK_RETURN:
            self._onLoseFocus(None)
        event.Skip()
                    
    def _onStringKey(self,event):
        event.Skip()
        if self.invalid: # check for validity after processing the keystroke
            wx.CallAfter(self.ShowStringValidity,True) # was invalid
        else:
            wx.CallAfter(self.ShowStringValidity,False) # was valid

    def _IndicateValidity(self):
        'Set the control colors to show invalid input'
        if self.invalid:
            self.SetForegroundColour("red")
            self.SetBackgroundColour("yellow")
            self.SetFocus()
            self.Refresh()
        else: # valid input
            self.SetBackgroundColour(
                wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
            self.SetForegroundColour("black")
            self.Refresh()

    def ShowStringValidity(self,previousInvalid=True):
        '''Check if input is valid. Anytime the input is
        invalid, call self.OKcontrol (if defined) because it is fast.
        If valid, check for any other invalid entries only when
        changing from invalid to valid, since that is slower.
        
        :param bool previousInvalid: True if the TextCtrl contents were
          invalid prior to the current change.
          
        '''
        val = self.GetValue().strip()
        self.invalid = not val
        self._IndicateValidity()
        if self.invalid:
            if self.OKcontrol:
                self.OKcontrol(False)
        elif self.OKcontrol and previousInvalid:
            self.OKcontrol(True)
        # always store the result
        if self.CIFinput: # for CIF make results ASCII
            self.result[self.key] = val.encode('ascii','replace') 
        else:
            self.result[self.key] = val
        log.LogVarChange(self.result,self.key)

    def _GetStringValue(self,event):
        '''Get string input and store.
        '''
        event.Skip() # process keystroke
        wx.CallAfter(self._SaveStringValue)
        
    def _SaveStringValue(self):
        val = self.GetValue().strip()
        # always store the result
        if self.CIFinput: # for CIF make results ASCII
            self.result[self.key] = val.encode('ascii','replace') 
        else:
            self.result[self.key] = val
        log.LogVarChange(self.result,self.key)

    def _onLoseFocus(self,event):
        if self.evaluated:
            self.EvaluateExpression()
        elif self.result is not None: # show formatted result, as Bob wants
            self._setValue(self.result[self.key])
        if self.OnLeave: self.OnLeave(invalid=self.invalid,
                                      value=self.result[self.key],
                                      tc=self,
                                      **self.OnLeaveArgs)
        if event: event.Skip()

    def EvaluateExpression(self):
        '''Show the computed value when an expression is entered to the TextCtrl
        Make sure that the number fits by truncating decimal places and switching
        to scientific notation, as needed. 
        Called on loss of focus, enter, etc..
        '''
        if self.invalid: return # don't substitute for an invalid expression
        if not self.evaluated: return # true when an expression is evaluated
        if self.result is not None: # retrieve the stored result
            self._setValue(self.result[self.key])
        self.evaluated = False # expression has been recast as value, reset flag
        
class NumberValidator(wx.PyValidator):
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

    :param number min: Minimum allowed value. If None (default) the
      lower limit is unbounded

    :param number max: Maximum allowed value. If None (default) the
      upper limit is unbounded
      
    :param dict/list result: List or dict where value should be placed when valid

    :param any key: key to use for result (int for list)

    :param function OKcontrol: function or class method to control
      an OK button for a window. 
      Ignored if None (default)

    :param bool CIFinput: allows use of a single '?' or '.' character
      as valid input.
      
    '''
    def __init__(self, typ, positiveonly=False, min=None, max=None,
                 result=None, key=None, OKcontrol=None, CIFinput=False):
        'Create the validator'
        wx.PyValidator.__init__(self)
        # save passed parameters
        self.typ = typ
        self.positiveonly = positiveonly
        self.min = min
        self.max = max
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
            self.validchars = '0123456789.-+eE/cosindcqrtap()*'
        else:
            self.validchars = None
            return
        if self.CIFinput:
            self.validchars += '?.'
    def Clone(self):
        'Create a copy of the validator, a strange, but required component'
        return NumberValidator(typ=self.typ, 
                               positiveonly=self.positiveonly,
                               min=self.min, max=self.max,
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
                log.LogVarChange(self.result,self.key)
                return
        try:
            val = self.typ(tc.GetValue())
        except (ValueError, SyntaxError) as e:
            if self.typ is float: # for float values, see if an expression can be evaluated
                val = G2py3.FormulaEval(tc.GetValue())
                if val is None:
                    tc.invalid = True
                    return
                else:
                    tc.evaluated = True
            else: 
                tc.invalid = True
                return
        # if self.max != None and self.typ == int:
        #     if val > self.max:
        #         tc.invalid = True
        # if self.min != None and self.typ == int:
        #     if val < self.min:
        #         tc.invalid = True  # invalid
        if self.max != None:
            if val > self.max:
                tc.invalid = True
        if self.min != None:
            if val < self.min:
                tc.invalid = True  # invalid
        if self.key is not None and self.result is not None and not tc.invalid:
            self.result[self.key] = val
            log.LogVarChange(self.result,self.key)

    def ShowValidity(self,tc):
        '''Set the control colors to show invalid input

        :param wx.TextCtrl tc: A reference to the TextCtrl that the validator
          is associated with.

        '''
        if tc.invalid:
            tc.SetForegroundColour("red")
            tc.SetBackgroundColour("yellow")
            tc.SetFocus()
            tc.Refresh()
            return False
        else: # valid input
            tc.SetBackgroundColour(
                wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
            tc.SetForegroundColour("black")
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
        if key == wx.WXK_RETURN:
            if tc.invalid:
                self.CheckInput(True) 
            else:
                self.CheckInput(False) 
            return
        if key < wx.WXK_SPACE or key == wx.WXK_DELETE or key > 255: # control characters get processed
            event.Skip()
            if tc.invalid:
                wx.CallAfter(self.CheckInput,True) 
            else:
                wx.CallAfter(self.CheckInput,False) 
            return
        elif chr(key) in self.validchars: # valid char pressed?
            event.Skip()
            if tc.invalid:
                wx.CallAfter(self.CheckInput,True) 
            else:
                wx.CallAfter(self.CheckInput,False) 
            return
        if not wx.Validator_IsSilent(): wx.Bell()
        return  # Returning without calling event.Skip, which eats the keystroke

class ASCIIValidator(wx.PyValidator):
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
        wx.PyValidator.__init__(self)
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
        self.result[self.key] = tc.GetValue().encode('ascii','replace')
        log.LogVarChange(self.result,self.key)

    def OnChar(self, event):
        '''Called each type a key is pressed
        ignores keys that are not allowed for int and float types
        '''
        key = event.GetKeyCode()
        tc = self.GetWindow()
        if key == wx.WXK_RETURN:
            self.TestValid(tc)
            return
        if key < wx.WXK_SPACE or key == wx.WXK_DELETE or key > 255: # control characters get processed
            event.Skip()
            self.TestValid(tc)
            return
        elif chr(key) in self.validchars: # valid char pressed?
            event.Skip()
            self.TestValid(tc)
            return
        if not wx.Validator_IsSilent():
            wx.Bell()
        return  # Returning without calling event.Skip, which eats the keystroke

################################################################################
#### Edit a large number of values
################################################################################
def CallScrolledMultiEditor(parent,dictlst,elemlst,prelbl=[],postlbl=[],
                 title='Edit items',header='',size=(300,250),
                             CopyButton=False, **kw):
    '''Shell routine to call a ScrolledMultiEditor dialog. See
    :class:`ScrolledMultiEditor` for parameter definitions.

    :returns: True if the OK button is pressed; False if the window is closed
      with the system menu or the Cancel button.

    '''
    dlg = ScrolledMultiEditor(parent,dictlst,elemlst,prelbl,postlbl,
                              title,header,size,
                              CopyButton, **kw)
    if dlg.ShowModal() == wx.ID_OK:
        dlg.Destroy()
        return True
    else:
        dlg.Destroy()
        return False

class ScrolledMultiEditor(wx.Dialog):
    '''Define a window for editing a potentially large number of dict- or
    list-contained values with validation for each item. Edited values are
    automatically placed in their source location. If invalid entries
    are provided, the TextCtrl is turned yellow and the OK button is disabled.

    The type for each TextCtrl validation is determined by the
    initial value of the entry (int, float or string). 
    Float values can be entered in the TextCtrl as numbers or also
    as algebraic expressions using operators + - / \* () and \*\*,
    in addition pi, sind(), cosd(), tand(), and sqrt() can be used,
    as well as appreviations s(), sin(), c(), cos(), t(), tan() and sq(). 

    :param wx.Frame parent: name of parent window, or may be None

    :param tuple dictlst: a list of dicts or lists containing values to edit

    :param tuple elemlst: a list of keys for each item in a dictlst. Must have the
      same length as dictlst.

    :param wx.Frame parent: name of parent window, or may be None
    
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
                 CopyButton=False,
                 minvals=[],maxvals=[],sizevals=[],
                 checkdictlst=[], checkelemlst=[], checklabel=""):
        if len(dictlst) != len(elemlst):
            raise Exception,"ScrolledMultiEditor error: len(dictlst) != len(elemlst) "+str(len(dictlst))+" != "+str(len(elemlst))
        if len(checkdictlst) != len(checkelemlst):
            raise Exception,"ScrolledMultiEditor error: len(checkdictlst) != len(checkelemlst) "+str(len(checkdictlst))+" != "+str(len(checkelemlst))
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
                if minvals[i] is not None: kargs['min']=minvals[i]
            if i < len(maxvals):
                if maxvals[i] is not None: kargs['max']=maxvals[i]
            if i < len(sizevals):
                if sizevals[i]: kargs['size']=sizevals[i]
            if CopyButton:
                import wx.lib.colourselect as wscs
                but = wscs.ColourSelect(label='v', # would like to use u'\u2193' or u'\u25BC' but not in WinXP
                                        # is there a way to test? 
                                        parent=panel,
                                        colour=(255,255,200),
                                        size=wx.Size(30,23),
                                        style=wx.RAISED_BORDER)
                but.Bind(wx.EVT_BUTTON, self._OnCopyButton)
                but.SetToolTipString('Press to copy adjacent value to all rows below')
                self.ButtonIndex[but] = i
                subSizer.Add(but)
            # create the validated TextCrtl, store it and add it to the sizer
            ctrl = ValidatedTxtCtrl(panel,d,k,OKcontrol=self.ControlOKButton,
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

################################################################################
#### 
################################################################################
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
    TreeItemType = G2frame.PatternTree.GetItemText(G2frame.PickId)
    copyopts = {'InTable':False,"startvalue":None,'currentsel':None}        
    hst = G2frame.PatternTree.GetItemText(G2frame.PatternId)
    histList = G2pdG.GetHistsLikeSelected(G2frame)
    if not histList:
        G2frame.ErrorDialog('No match','No histograms match '+hst,G2frame.dataFrame)
        return
    dlg = wx.Dialog(G2frame.dataDisplay,wx.ID_ANY,'Set a parameter value',
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
    dlg = G2MultiChoiceDialog(
        G2frame.dataFrame, 
        'Copy parameter '+lbl+' from\n'+hst,
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
    Id = GetPatternTreeItemId(G2frame,G2frame.root,hst)
    hstData = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'))[0]
    for h in copyList:
        Id = GetPatternTreeItemId(G2frame,G2frame.root,h)
        instData = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,Id,'Instrument Parameters'))[0]
        if len(hstData) != len(instData) or hstData['Type'][0] != instData['Type'][0]:  #don't mix data types or lam & lam1/lam2 parms!
            print h+' not copied - instrument parameters not commensurate'
            continue
        hData = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,Id,TreeItemType))
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
            G2frame.dataDisplay,dictlst,
            len(dictlst)*[keyLst[-1]],prelbl,
            header='Editing parameter '+lbl,
            CopyButton=True,**args)
        dlg.CenterOnParent()
        if dlg.ShowModal() != wx.ID_OK:
            array.update(saveArray)
        dlg.Destroy()

################################################################################
#### Custom checkbox that saves values into dict/list as used
################################################################################
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
    '''
    def __init__(self,parent,label,loc,key):
        wx.CheckBox.__init__(self,parent,id=wx.ID_ANY,label=label)
        self.loc = loc
        self.key = key
        self.SetValue(self.loc[self.key]==True)
        self.Bind(wx.EVT_CHECKBOX, self._OnCheckBox)
    def _OnCheckBox(self,event):
        self.loc[self.key] = self.GetValue()
        log.LogVarChange(self.loc,self.key)

################################################################################
####
################################################################################
class PickTwoDialog(wx.Dialog):
    '''This does not seem to be in use
    '''
    def __init__(self,parent,title,prompt,names,choices):
        wx.Dialog.__init__(self,parent,-1,title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
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
            
        self.panel.DestroyChildren()
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
        parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

################################################################################
#### Column-order selection
################################################################################

def GetItemOrder(parent,keylist,vallookup,posdict):
    '''Creates a panel where items can be ordered into columns
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
    sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.EXPAND|wx.ALL, 5)
    dlg.SetSizer(sizer)
    sizer.Fit(dlg)
    val = dlg.ShowModal()

################################################################################
####
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
        self.SetBackgroundColour(WHITE)
        self.SetSizer(self.GBsizer)
        colList = [str(i) for i in range(self.maxcol+2)]
        for i in range(self.maxcol+1):
            wid = wx.StaticText(self,wx.ID_ANY,str(i),style=wx.ALIGN_CENTER)
            wid.SetBackgroundColour(DULL_YELLOW)
            wid.SetMinSize((50,-1))
            self.GBsizer.Add(wid,(0,i),flag=wx.ALIGN_CENTER|wx.EXPAND)
        self.chceDict = {}
        for row,nam in enumerate(self.keylist):
            posdict = self.posdict[nam]
            for col in posdict:
                lbl = posdict[col]
                pnl = wx.Panel(self,wx.ID_ANY)
                pnl.SetBackgroundColour(VERY_LIGHT_GREY)
                insize = wx.BoxSizer(wx.VERTICAL)
                wid = wx.Choice(pnl,wx.ID_ANY,choices=colList)
                insize.Add(wid,0,wx.EXPAND|wx.BOTTOM,3)
                wid.SetSelection(col)
                self.chceDict[wid] = (row,col)
                wid.Bind(wx.EVT_CHOICE,self.OnChoice)
                wid = wx.StaticText(pnl,wx.ID_ANY,lbl)
                insize.Add(wid,0,flag=wx.EXPAND)
                val = G2py3.FormatSigFigs(self.vallookup[nam][lbl],maxdigits=8)
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
                self.GBsizer.Add(wid,(0,i),flag=wx.ALIGN_CENTER|wx.EXPAND)
            colList = [str(i) for i in range(self.maxcol+2)]
            for wid in self.chceDict:
                wid.SetItems(colList)
                wid.SetSelection(self.chceDict[wid][1])
        self.GBsizer.Layout()
        self.FitInside()

################################################################################
#### Help support routines
################################################################################
################################################################################
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
    def __init__(self,frame,helpType=None,helpLbl=None,morehelpitems=[],title=''):
        wx.Menu.__init__(self,title)
        self.HelpById = {}
        self.frame = frame
        self.Append(help='', id=wx.ID_ABOUT, kind=wx.ITEM_NORMAL,
            text='&About GSAS-II')
        frame.Bind(wx.EVT_MENU, self.OnHelpAbout, id=wx.ID_ABOUT)
        if GSASIIpath.whichsvn():
            helpobj = self.Append(
                help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
                text='&Check for updates')
            frame.Bind(wx.EVT_MENU, self.OnCheckUpdates, helpobj)
            helpobj = self.Append(
                help='', id=wx.ID_ANY, kind=wx.ITEM_NORMAL,
                text='&Regress to an old GSAS-II version')
            frame.Bind(wx.EVT_MENU, self.OnSelectVersion, helpobj)
        for lbl,indx in morehelpitems:
            helpobj = self.Append(text=lbl,
                id=wx.ID_ANY, kind=wx.ITEM_NORMAL)
            frame.Bind(wx.EVT_MENU, self.OnHelpById, helpobj)
            self.HelpById[helpobj.GetId()] = indx
        # add a help item only when helpType is specified
        if helpType is not None:
            self.AppendSeparator()
            if helpLbl is None: helpLbl = helpType
            helpobj = self.Append(text='Help on '+helpLbl,
                                  id=wx.ID_ANY, kind=wx.ITEM_NORMAL)
            frame.Bind(wx.EVT_MENU, self.OnHelpById, helpobj)
            self.HelpById[helpobj.GetId()] = helpType
       
    def OnHelpById(self,event):
        '''Called when Help on... is pressed in a menu. Brings up
        a web page for documentation.
        '''
        helpType = self.HelpById.get(event.GetId())
        if helpType is None:
            print 'Error: help lookup failed!',event.GetEventObject()
            print 'id=',event.GetId()
        elif helpType == 'OldTutorials': # this will go away
            self.frame.Tutorials = True
            ShowHelp(helpType,self.frame)
        elif helpType == 'Tutorials': 
            dlg = OpenTutorial(self.frame)
            if dlg.ShowModal() == wx.ID_OK:
                self.frame.Tutorials = True
            dlg.Destroy()
            return
        else:
            ShowHelp(helpType,self.frame)

    def OnHelpAbout(self, event):
        "Display an 'About GSAS-II' box"
        import GSASII
        info = wx.AboutDialogInfo()
        info.Name = 'GSAS-II'
        ver = GSASIIpath.svnGetRev()
        if ver: 
            info.Version = 'Revision '+str(ver)+' (svn), version '+GSASII.__version__
        else:
            info.Version = 'Revision '+str(GSASIIpath.GetVersionNumber())+' (.py files), version '+GSASII.__version__
        #info.Developers = ['Robert B. Von Dreele','Brian H. Toby']
        info.Copyright = ('(c) ' + time.strftime('%Y') +
''' Argonne National Laboratory
This product includes software developed
by the UChicago Argonne, LLC, as 
Operator of Argonne National Laboratory.''')
        info.Description = '''General Structure Analysis System-II (GSAS-II)
Robert B. Von Dreele and Brian H. Toby

Please cite as:
B.H. Toby & R.B. Von Dreele, J. Appl. Cryst. 46, 544-549 (2013) '''

        info.WebSite = ("https://subversion.xray.aps.anl.gov/trac/pyGSAS","GSAS-II home page")
        wx.AboutBox(info)

    def OnCheckUpdates(self,event):
        '''Check if the GSAS-II repository has an update for the current source files
        and perform that update if requested.
        '''
        if not GSASIIpath.whichsvn():
            dlg = wx.MessageDialog(self.frame,
                                   'No Subversion','Cannot update GSAS-II because subversion (svn) was not found.',
                                   wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return
        wx.BeginBusyCursor()
        local = GSASIIpath.svnGetRev()
        if local is None: 
            wx.EndBusyCursor()
            dlg = wx.MessageDialog(self.frame,
                                   'Unable to run subversion on the GSAS-II current directory. Is GSAS-II installed correctly?',
                                   'Subversion error',
                                   wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return
        print 'Installed GSAS-II version: '+local
        repos = GSASIIpath.svnGetRev(local=False)
        wx.EndBusyCursor()
        if repos is None: 
            dlg = wx.MessageDialog(self.frame,
                                   'Unable to access the GSAS-II server. Is this computer on the internet?',
                                   'Server unavailable',
                                   wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return
        print 'GSAS-II version on server: '+repos
        if local == repos:
            dlg = wx.MessageDialog(self.frame,
                                   'GSAS-II is up-to-date. Version '+local+' is already loaded.',
                                   'GSAS-II Up-to-date',
                                   wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return
        mods = GSASIIpath.svnFindLocalChanges()
        if mods:
            dlg = wx.MessageDialog(self.frame,
                                   'You have version '+local+
                                   ' of GSAS-II installed, but the current version is '+repos+
                                   '. However, '+str(len(mods))+
                                   ' file(s) on your local computer have been modified.'
                                   ' Updating will attempt to merge your local changes with '
                                   'the latest GSAS-II version, but if '
                                   'conflicts arise, local changes will be '
                                   'discarded. It is also possible that the '
                                   'local changes my prevent GSAS-II from running. '
                                   'Press OK to start an update if this is acceptable:',
                                   'Local GSAS-II Mods',
                                   wx.OK|wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return
            else:
                dlg.Destroy()
        else:
            dlg = wx.MessageDialog(self.frame,
                                   'You have version '+local+
                                   ' of GSAS-II installed, but the current version is '+repos+
                                   '. Press OK to start an update:',
                                   'GSAS-II Updates',
                                   wx.OK|wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return
            dlg.Destroy()
        print 'start updates'
        dlg = wx.MessageDialog(self.frame,
                               'Your project will now be saved, GSAS-II will exit and an update '
                               'will be performed and GSAS-II will restart. Press Cancel to '
                               'abort the update',
                               'Start update?',
                               wx.OK|wx.CANCEL)
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        dlg.Destroy()
        self.frame.OnFileSave(event)
        GSASIIpath.svnUpdateProcess(projectfile=self.frame.GSASprojectfile)
        return

    def OnSelectVersion(self,event):
        '''Allow the user to select a specific version of GSAS-II
        '''
        if not GSASIIpath.whichsvn():
            dlg = wx.MessageDialog(self,'No Subversion','Cannot update GSAS-II because subversion (svn) '+
                                   'was not found.'
                                   ,wx.OK)
            dlg.ShowModal()
            return
        local = GSASIIpath.svnGetRev()
        if local is None: 
            dlg = wx.MessageDialog(self.frame,
                                   'Unable to run subversion on the GSAS-II current directory. Is GSAS-II installed correctly?',
                                   'Subversion error',
                                   wx.OK)
            dlg.ShowModal()
            return
        mods = GSASIIpath.svnFindLocalChanges()
        if mods:
            dlg = wx.MessageDialog(self.frame,
                                   'You have version '+local+
                                   ' of GSAS-II installed'
                                   '. However, '+str(len(mods))+
                                   ' file(s) on your local computer have been modified.'
                                   ' Downdating will attempt to merge your local changes with '
                                   'the selected GSAS-II version. '
                                   'Downdating is not encouraged because '
                                   'if merging is not possible, your local changes will be '
                                   'discarded. It is also possible that the '
                                   'local changes my prevent GSAS-II from running. '
                                   'Press OK to continue anyway.',
                                   'Local GSAS-II Mods',
                                   wx.OK|wx.CANCEL)
            if dlg.ShowModal() != wx.ID_OK:
                dlg.Destroy()
                return
            dlg.Destroy()
        dlg = downdate(parent=self.frame)
        if dlg.ShowModal() == wx.ID_OK:
            ver = dlg.getVersion()
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        print('start regress to '+str(ver))
        GSASIIpath.svnUpdateProcess(
            projectfile=self.frame.GSASprojectfile,
            version=str(ver)
            )
        self.frame.OnFileSave(event)
        return

################################################################################
class AddHelp(wx.Menu):
    '''For the Mac: creates an entry to the help menu of type 
    'Help on <helpType>': where helpType is a reference to an HTML page to
    be opened.

    NOTE: when appending this menu (menu.Append) be sure to set the title to
    '&Help' so that wx handles it correctly.
    '''
    def __init__(self,frame,helpType,helpLbl=None,title=''):
        wx.Menu.__init__(self,title)
        self.frame = frame
        if helpLbl is None: helpLbl = helpType
        # add a help item only when helpType is specified
        helpobj = self.Append(text='Help on '+helpLbl,
                              id=wx.ID_ANY, kind=wx.ITEM_NORMAL)
        frame.Bind(wx.EVT_MENU, self.OnHelpById, helpobj)
        self.HelpById = helpType
       
    def OnHelpById(self,event):
        '''Called when Help on... is pressed in a menu. Brings up
        a web page for documentation.
        '''
        ShowHelp(self.HelpById,self.frame)

################################################################################
class HelpButton(wx.Button):
    '''Create a help button that displays help information.
    The text is displayed in a modal message window.

    TODO: it might be nice if it were non-modal: e.g. it stays around until
    the parent is deleted or the user closes it, but this did not work for
    me. 

    :param parent: the panel which will be the parent of the button
    :param str msg: the help text to be displayed
    '''
    def __init__(self,parent,msg):
        if sys.platform == "darwin": 
            wx.Button.__init__(self,parent,wx.ID_HELP)
        else:
            wx.Button.__init__(self,parent,wx.ID_ANY,'?',style=wx.BU_EXACTFIT)
        self.Bind(wx.EVT_BUTTON,self._onPress)
        self.msg=StripIndents(msg)
        self.parent = parent
    def _onClose(self,event):
        self.dlg.EndModal(wx.ID_CANCEL)
    def _onPress(self,event):
        'Respond to a button press by displaying the requested text'
        #dlg = wx.MessageDialog(self.parent,self.msg,'Help info',wx.OK)
        self.dlg = wx.Dialog(self.parent,wx.ID_ANY,'Help information', 
                        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        #self.dlg.SetBackgroundColour(wx.WHITE)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        txt = wx.StaticText(self.dlg,wx.ID_ANY,self.msg)
        mainSizer.Add(txt,1,wx.ALL|wx.EXPAND,10)
        txt.SetBackgroundColour(wx.WHITE)

        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(self.dlg, wx.ID_CLOSE) 
        btn.Bind(wx.EVT_BUTTON,self._onClose)
        btnsizer.Add(btn)
        mainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        self.dlg.SetSizer(mainSizer)
        mainSizer.Fit(self.dlg)
        self.dlg.ShowModal()
        self.dlg.Destroy()
################################################################################
class MyHtmlPanel(wx.Panel):
    '''Defines a panel to display HTML help information, as an alternative to
    displaying help information in a web browser.
    '''
    def __init__(self, frame, id):
        self.frame = frame
        wx.Panel.__init__(self, frame, id)
        sizer = wx.BoxSizer(wx.VERTICAL)
        back = wx.Button(self, -1, "Back")
        back.Bind(wx.EVT_BUTTON, self.OnBack)
        self.htmlwin = G2HtmlWindow(self, id, size=(750,450))
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
            event.Skip()
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
def StripIndents(msg):
    'Strip indentation from multiline strings'
    msg1 = msg.replace('\n ','\n')
    while msg != msg1:
        msg = msg1
        msg1 = msg.replace('\n ','\n')
    return msg.replace('\n\t','\n')

def G2MessageBox(parent,msg,title='Error'):
    '''Simple code to display a error or warning message
    '''
    dlg = wx.MessageDialog(parent,StripIndents(msg), title, wx.OK)
    dlg.ShowModal()
    dlg.Destroy()
        
################################################################################
class downdate(wx.Dialog):
    '''Dialog to allow a user to select a version of GSAS-II to install
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
        sizer.Add(line, 0, wx.EXPAND|wx.ALIGN_CENTER|wx.ALL, 10)

        self.text = wx.StaticText(pnl,  wx.ID_ANY, "")
        sizer.Add(self.text, 0, wx.ALIGN_LEFT|wx.EXPAND|wx.ALL, 5)

        line = wx.StaticLine(pnl,-1, size=(-1,3), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.EXPAND|wx.ALIGN_CENTER|wx.ALL, 10)
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
#### Display Help information
################################################################################
# define some globals 
htmlPanel = None
htmlFrame = None
htmlFirstUse = True
helpLocDict = {}
path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # save location of this file
def ShowHelp(helpType,frame):
    '''Called to bring up a web page for documentation.'''
    global htmlFirstUse
    # look up a definition for help info from dict
    helplink = helpLocDict.get(helpType)
    if helplink is None:
        # no defined link to use, create a default based on key
        helplink = 'gsasII.html#'+helpType.replace(' ','_')
    helplink = os.path.join(path2GSAS2,'help',helplink)
    # determine if a web browser or the internal viewer should be used for help info
    if GSASIIpath.GetConfigValue('Help_mode'):
        helpMode = GSASIIpath.GetConfigValue('Help_mode')
    else:
        helpMode = 'browser'
    if helpMode == 'internal':
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
        pfx = "file://"
        if sys.platform.lower().startswith('win'):
            pfx = ''
        if htmlFirstUse:
            webbrowser.open_new(pfx+helplink)
            htmlFirstUse = False
        else:
            webbrowser.open(pfx+helplink, new=0, autoraise=True)
def ShowWebPage(URL,frame):
    '''Called to show a tutorial web page.
    '''
    global htmlFirstUse
    # determine if a web browser or the internal viewer should be used for help info
    if GSASIIpath.GetConfigValue('Help_mode'):
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
################################################################################
#### Tutorials selector
################################################################################
G2BaseURL = "https://subversion.xray.aps.anl.gov/pyGSAS"
# N.B. tutorialCatalog is generated by routine catalog.py, which also generates the appropriate
# empty directories (.../MT/* .../trunk/GSASII/* *=[help,Exercises])
tutorialCatalog = (
    # tutorial dir,      exercise dir,      web page file name                                title for page

    ['StartingGSASII', 'StartingGSASII', 'Starting GSAS.htm',
       'Starting GSAS-II'],
       
    ['FitPeaks', 'FitPeaks', 'Fit Peaks.htm',
       'Fitting individual peaks & autoindexing'],
       
    ['CWNeutron', 'CWNeutron', 'Neutron CW Powder Data.htm',
       'CW Neutron Powder fit for Yttrium-Iron Garnet'],
    ['LabData', 'LabData', 'Laboratory X.htm',
       'Fitting laboratory X-ray powder data for fluoroapatite'],
    ['CWCombined', 'CWCombined', 'Combined refinement.htm',
       'Combined X-ray/CW-neutron refinement of PbSO4'],
    ['TOF-CW Joint Refinement', 'TOF-CW Joint Refinement', 'TOF combined XN Rietveld refinement in GSAS.htm',
       'Combined X-ray/TOF-neutron Rietveld refinement'],
    ['SeqRefine', 'SeqRefine', 'SequentialTutorial.htm',
       'Sequential refinement of multiple datasets'],
    ['SeqParametric', 'SeqParametric', 'ParametricFitting.htm',
       'Parametric Fitting and Pseudo Variables for Sequential Fits'],
       
    ['CFjadarite', 'CFjadarite', 'Charge Flipping in GSAS.htm',
       'Charge Flipping structure solution for jadarite'],
    ['CFsucrose', 'CFsucrose', 'Charge Flipping - sucrose.htm',
       'Charge Flipping structure solution for sucrose'],
    ['TOF Charge Flipping', 'TOF Charge Flipping', 'Charge Flipping with TOF single crystal data in GSASII.htm',
       'Charge flipping with neutron TOF single crystal data'],
    ['MCsimanneal', 'MCsimanneal', 'MCSA in GSAS.htm',
       'Monte-Carlo simulated annealing structure'],

    ['2DCalibration', '2DCalibration', 'Calibration of an area detector in GSAS.htm',
       'Calibration of an area detector'],
    ['2DIntegration', '2DIntegration', 'Integration of area detector data in GSAS.htm',
       'Integration of area detector data'],
    ['TOF Calibration', 'TOF Calibration', 'Calibration of a TOF powder diffractometer.htm',
       'Calibration of a Neutron TOF diffractometer'],
       
    ['2DStrain', '2DStrain', 'Strain fitting of 2D data in GSAS-II.htm',
       'Strain fitting of 2D data'],
       
    ['SAimages', 'SAimages', 'Small Angle Image Processing.htm',
       'Image Processing of small angle x-ray data'],
    ['SAfit', 'SAfit', 'Fitting Small Angle Scattering Data.htm',
       'Fitting small angle x-ray data (alumina powder)'],
    ['SAsize', 'SAsize', 'Small Angle Size Distribution.htm',
       'Small angle x-ray data size distribution (alumina powder)'],
    ['SAseqref', 'SAseqref', 'Sequential Refinement of Small Angle Scattering Data.htm',
       'Sequential refinement with small angle scattering data'],
    
    #['TOF Sequential Single Peak Fit', 'TOF Sequential Single Peak Fit', '', ''],
    #['TOF Single Crystal Refinement', 'TOF Single Crystal Refinement', '', ''],
    )
if GSASIIpath.GetConfigValue('Tutorial_location'):
    tutorialPath = GSASIIpath.GetConfigValue('Tutorial_location')
else:
    tutorialPath = GSASIIpath.path2GSAS2

class OpenTutorial(wx.Dialog):
    '''Open a tutorial, optionally copying it to the local disk. Always copy
    the data files locally.

    For now tutorials will always be copied into the source code tree, but it
    might be better to have an option to copy them somewhere else, for people
    who don't have write access to the GSAS-II source code location. 
    '''
    # TODO: set default input-file open location to the download location
    def __init__(self,parent=None):
        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, wx.ID_ANY, 'Open Tutorial', style=style)
        self.frame = parent
        pnl = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(
            pnl,  wx.ID_ANY,
            'Select the tutorial to be run and the mode of access'
            )
        sizer.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        msg = '''To save download time for GSAS-II tutorials and their
        sample data files are being moved out of the standard
        distribution. This dialog allows users to load selected
        tutorials to their computer.

        Tutorials can be viewed over the internet or downloaded
        to this computer. The sample data can be downloaded or not,
        (but it is not possible to run the tutorial without the
        data). If no web access is available, tutorials that were
        previously downloaded can be viewed.

        By default, files are downloaded into the location used
        for the GSAS-II distribution, but this may not be possible
        if the software is installed by a administrator. The
        download location can be changed using the "Set data
        location" or the "Tutorial_location" configuration option
        (see config_example.py).
        '''
        hlp = HelpButton(pnl,msg)
        sizer.Add(hlp,0,wx.ALIGN_RIGHT|wx.ALL)
        #======================================================================
        # # This is needed only until we get all the tutorials items moved
        # btn = wx.Button(pnl, wx.ID_ANY, "Open older tutorials") 
        # btn.Bind(wx.EVT_BUTTON, self.OpenOld)
        # sizer.Add(btn,0,wx.ALIGN_CENTRE|wx.ALL)
        #======================================================================
        self.BrowseMode = 1
        choices = [
            'make local copy of tutorial and data, then open',
            'run from web (copy data locally)',
            'browse on web (data not loaded)', 
            'open from local tutorial copy',
        ]
        self.mode = wx.RadioBox(pnl,wx.ID_ANY,'access mode:',
                                wx.DefaultPosition, wx.DefaultSize,
                                choices, 1, wx.RA_SPECIFY_COLS)
        self.mode.SetSelection(self.BrowseMode)
        self.mode.Bind(wx.EVT_RADIOBOX, self.OnModeSelect)
        sizer.Add(self.mode,0,WACV)
        sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(pnl, wx.ID_ANY, "Set download location")
        btn.Bind(wx.EVT_BUTTON, self.SelectDownloadLoc)
        sizer1.Add(btn,0,WACV)
        self.dataLoc = wx.StaticText(pnl, wx.ID_ANY,tutorialPath)
        sizer1.Add(self.dataLoc,0,WACV)
        sizer.Add(sizer1)
        self.listbox = wx.ListBox(pnl, wx.ID_ANY, size=(450, 100), style=wx.LB_SINGLE)
        self.listbox.Bind(wx.EVT_LISTBOX, self.OnTutorialSelected)
        self.OnModeSelect(None)
        #self.FillListBox()
        sizer.Add(self.listbox,1,WACV|wx.EXPAND|wx.ALL,1)
        
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(pnl, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        pnl.SetSizer(sizer)
        sizer.Fit(self)
        self.topsizer=sizer
        self.CenterOnParent()
    # def OpenOld(self,event):
    #     '''Open old tutorials. This is needed only until we get all the tutorials items moved
    #     '''
    #     self.EndModal(wx.ID_OK)
    #     self.frame.Tutorials = True 
    #     ShowHelp('Tutorials',self.frame)
    def OnModeSelect(self,event):
        '''Respond when the mode is changed
        '''
        self.BrowseMode = self.mode.GetSelection()
        if self.BrowseMode == 3:
            import glob
            filelist = glob.glob(os.path.join(tutorialPath,'help','*','*.htm'))
            taillist = [os.path.split(f)[1] for f in filelist]
            itemlist = [tut[-1] for tut in tutorialCatalog if tut[2] in taillist]
        else:
            itemlist = [tut[-1] for tut in tutorialCatalog if tut[-1]]
        self.listbox.Clear()
        self.listbox.AppendItems(itemlist)
    def OnTutorialSelected(self,event):
        '''Respond when a tutorial is selected. Load tutorials and data locally,
        as needed and then display the page
        '''
        for tutdir,exedir,htmlname,title in tutorialCatalog:
            if title == event.GetString(): break
        else:
            raise Exception("Match to file not found")
        if self.BrowseMode == 0 or self.BrowseMode == 1:
            try: 
                self.ValidateTutorialDir(tutorialPath,G2BaseURL)
            except:
                G2MessageBox(self.frame,
            '''The selected directory is not valid.
            
            You must use a directory that you have write access
            to. You can reuse a directory previously used for 
            downloads, but the help and Tutorials subdirectories
             must be created by this routine. 
            ''')
                return
        self.dataLoc.SetLabel(tutorialPath)
        self.EndModal(wx.ID_OK)
        wx.BeginBusyCursor()
        if self.BrowseMode == 0:
            # xfer data & web page locally, then open web page
            self.LoadTutorial(tutdir,tutorialPath,G2BaseURL)
            self.LoadExercise(exedir,tutorialPath,G2BaseURL)
            URL = os.path.join(tutorialPath,'help',tutdir,htmlname)
            ShowWebPage(URL,self.frame)
        elif self.BrowseMode == 1:
            # xfer data locally, open web page remotely
            self.LoadExercise(exedir,tutorialPath,G2BaseURL)
            URL = os.path.join(G2BaseURL,'Tutorials',tutdir,htmlname)
            ShowWebPage(URL,self.frame)
        elif self.BrowseMode == 2:
            # open web page remotely, don't worry about data
            URL = os.path.join(G2BaseURL,'Tutorials',tutdir,htmlname)
            ShowWebPage(URL,self.frame)
        elif self.BrowseMode == 3:
            # open web page that has already been transferred
            URL = os.path.join(tutorialPath,'help',tutdir,htmlname)
            ShowWebPage(URL,self.frame)
        else:
            wx.EndBusyCursor()
            raise Exception("How did this happen!")
        wx.EndBusyCursor()
    def ValidateTutorialDir(self,fullpath=tutorialPath,baseURL=G2BaseURL):
        '''Load help to new directory or make sure existing directory looks correctly set up
        throws an exception if there is a problem.
        '''
        if os.path.exists(fullpath):
            if os.path.exists(os.path.join(fullpath,"help")):
                if not GSASIIpath.svnGetRev(os.path.join(fullpath,"help")):
                    print("Problem with "+fullpath+" dir help exists but is not in SVN")
                    raise Exception
            if os.path.exists(os.path.join(fullpath,"Exercises")):
                if not GSASIIpath.svnGetRev(os.path.join(fullpath,"Exercises")):
                    print("Problem with "+fullpath+" dir Exercises exists but is not in SVN")
                    raise Exception
            if (os.path.exists(os.path.join(fullpath,"help")) and
                    os.path.exists(os.path.join(fullpath,"Exercises"))):
                return True # both good
            elif (os.path.exists(os.path.join(fullpath,"help")) or
                    os.path.exists(os.path.join(fullpath,"Exercises"))):
                print("Problem: dir "+fullpath+" exists has either help or Exercises, not both")
                raise Exception
        wx.BeginBusyCursor()
        if not GSASIIpath.svnInstallDir(baseURL+"/MT",fullpath):
            wx.EndBusyCursor()
            print("Problem transferring empty directory from web")
            raise Exception
        wx.EndBusyCursor()
        return True

    def LoadTutorial(self,tutorialname,fullpath=tutorialPath,baseURL=G2BaseURL):
        'Load a Tutorial to the selected location'
        if GSASIIpath.svnSwitchDir("help",tutorialname,baseURL+"/Tutorials",fullpath):
            return True
        print("Problem transferring Tutorial from web")
        raise Exception
        
    def LoadExercise(self,tutorialname,fullpath=tutorialPath,baseURL=G2BaseURL):
        'Load Exercise file(s) for a Tutorial to the selected location'
        if GSASIIpath.svnSwitchDir("Exercises",tutorialname,baseURL+"/Exercises",fullpath):
            return True
        print ("Problem transferring Exercise from web")
        raise Exception
        
    def SelectDownloadLoc(self,event):
        '''Select a download location,
        Cancel resets to the default
        '''
        global tutorialPath
        localpath = os.path.abspath(os.path.expanduser('~/G2tutorials'))
        dlg = wx.DirDialog(self, "Choose a directory for downloads:",
                           defaultPath=localpath)#,style=wx.DD_DEFAULT_STYLE)
                           #)
        if dlg.ShowModal() == wx.ID_OK:
            pth = dlg.GetPath()
        else:
            if GSASIIpath.GetConfigValue('Tutorial_location'):
                pth = GSASIIpath.GetConfigValue('Tutorial_location')
            else:
                pth = GSASIIpath.path2GSAS2
        if not os.path.exists(pth):
            try:
                os.makedirs(pth)
            except OSError:
                msg = 'The selected directory is not valid.\n\t'
                msg += pth
                msg += '\n\nAn attempt to create the directory failed'
                G2MessageBox(self.frame,msg)
                return
        self.ValidateTutorialDir(pth,G2BaseURL)
        try: 
            self.ValidateTutorialDir(pth,G2BaseURL)
            tutorialPath = pth
        except:
            G2MessageBox(self.frame,
            '''The selected directory is not valid.
            
            You must use a directory that you have write access
            to. You can reuse a directory previously used for 
            downloads, but the help and Tutorials subdirectories
             must be created by this routine. 
            ''')
        self.dataLoc.SetLabel(tutorialPath)
    
if __name__ == '__main__':
    app = wx.PySimpleApp()
    GSASIIpath.InvokeDebugOpts()
    frm = wx.Frame(None) # create a frame
    frm.Show(True)
    dlg = OpenTutorial(frm)
    if dlg.ShowModal() == wx.ID_OK:
        print "OK"
    else:
        print "Cancel"
    dlg.Destroy()
    import sys
    sys.exit()
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
    Data3 = {
         'Order':1.0,
         'omega':1.1,
         'chi':2.0,
         'phi':2.3,
         'Order1':1.0,
         'omega1':1.1,
         'chi1':2.0,
         'phi1':2.3,
         'Order2':1.0,
         'omega2':1.1,
         'chi2':2.0,
         'phi2':2.3,
         }
    elemlst = sorted(Data3.keys())
    dictlst = len(elemlst)*[Data3,]
    prelbl = elemlst[:]
    prelbl[0]="this is a much longer label to stretch things out"
    Data2 = len(elemlst)*[False,]
    Data2[1] = Data2[3] = True
    Checkdictlst = len(elemlst)*[Data2,]
    Checkelemlst = range(len(Checkdictlst))
    #print 'before',Data3,'\n',Data2
    #print dictlst,"\n",elemlst
    #print Checkdictlst,"\n",Checkelemlst
    dlg = ScrolledMultiEditor(
        frm,dictlst,elemlst,prelbl,
        checkdictlst=Checkdictlst,checkelemlst=Checkelemlst,
        checklabel="Refine?",
        header="test",CopyButton=True)
    if dlg.ShowModal() == wx.ID_OK:
        print "OK"
    else:
        print "Cancel"
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
    # dlg = G2MultiChoiceDialog(frm, 'Sequential refinement',
    #                           'Select dataset to include',
    #                           choices)
    # sel = range(2,11,2)
    # dlg.SetSelections(sel)
    # dlg.SetSelections((1,5))
    # if dlg.ShowModal() == wx.ID_OK:
    #     for sel in dlg.GetSelections():
    #         print sel,choices[sel]
    
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

    pnl = wx.Panel(frm)
    siz = wx.BoxSizer(wx.VERTICAL)

    td = {'Goni':200.,'a':1.,'calc':1./3.,'string':'s'}
    for key in sorted(td):
        txt = ValidatedTxtCtrl(pnl,td,key)
        siz.Add(txt)
    pnl.SetSizer(siz)
    siz.Fit(frm)
    app.MainLoop()
    print td
