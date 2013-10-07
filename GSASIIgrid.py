# -*- coding: utf-8 -*-
#GSASIIgrid - data display routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIgrid: Basic GUI routines*
--------------------------------

'''
import wx
import wx.grid as wg
import wx.wizard as wz
import wx.aui
import wx.lib.scrolledpanel as wxscroll
import time
import cPickle
import sys
import numpy as np
import numpy.ma as ma
import os.path
import wx.html        # could postpone this for quicker startup
import webbrowser     # could postpone this for quicker startup
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIImath as G2mth
import GSASIIIO as G2IO
import GSASIIlattice as G2lat
import GSASIIplot as G2plt
import GSASIIpwdGUI as G2pdG
import GSASIIimgGUI as G2imG
import GSASIIphsGUI as G2phG
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIIconstrGUI as G2cnstG
import GSASIIrestrGUI as G2restG
import GSASIIpy3 as G2py3

# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)

# globals we will use later
__version__ = None # gets overridden in GSASII.py
path2GSAS2 = os.path.dirname(os.path.realpath(__file__)) # save location of this file
helpLocDict = {}
htmlPanel = None
htmlFrame = None
helpMode = 'browser'
#if sys.platform.lower().startswith('win'): helpMode = 'internal' # need a global control to set this
    
htmlFirstUse = True

[ wxID_FOURCALC, wxID_FOURSEARCH, wxID_FOURCLEAR, wxID_PEAKSMOVE, wxID_PEAKSCLEAR, 
    wxID_CHARGEFLIP, wxID_PEAKSUNIQUE, wxID_PEAKSDELETE, wxID_PEAKSDA,
    wxID_PEAKSDISTVP, wxID_PEAKSVIEWPT, wxID_FINDEQVPEAKS,wxID_SHOWBONDS,wxID_MULTIMCSA,
    wxID_SINGLEMCSA
] = [wx.NewId() for item in range(15)]

[ wxID_PWDRADD, wxID_HKLFADD,wxID_PWDANALYSIS,wxID_DATADELETE,
] = [wx.NewId() for item in range(4)]

[ wxID_ATOMSEDITADD, wxID_ATOMSEDITINSERT, wxID_ATOMSEDITDELETE, wxID_ATOMSREFINE, 
    wxID_ATOMSMODIFY, wxID_ATOMSTRANSFORM, wxID_ATOMSVIEWADD, wxID_ATOMVIEWINSERT,
    wxID_RELOADDRAWATOMS,wxID_ATOMSDISAGL,wxID_ATOMMOVE,
    wxID_ASSIGNATMS2RB,wxID_ATOMSPDISAGL
] = [wx.NewId() for item in range(13)]

[ wxID_DRAWATOMSTYLE, wxID_DRAWATOMLABEL, wxID_DRAWATOMCOLOR, wxID_DRAWATOMRESETCOLOR, 
    wxID_DRAWVIEWPOINT, wxID_DRAWTRANSFORM, wxID_DRAWDELETE, wxID_DRAWFILLCELL, 
    wxID_DRAWADDEQUIV, wxID_DRAWFILLCOORD, wxID_DRAWDISAGLTOR,  wxID_DRAWPLANE,
    wxID_DRAWDISTVP,
] = [wx.NewId() for item in range(13)]

[ wxID_DRAWRESTRBOND, wxID_DRAWRESTRANGLE, wxID_DRAWRESTRPLANE, wxID_DRAWRESTRCHIRAL,
] = [wx.NewId() for item in range(4)]

[ wxID_ADDMCSAATOM,wxID_ADDMCSARB,wxID_CLEARMCSARB,wxID_MOVEMCSA,wxID_MCSACLEARRESULTS,
] = [wx.NewId() for item in range(5)]

[ wxID_CLEARTEXTURE,wxID_REFINETEXTURE,
] = [wx.NewId() for item in range(2)]

[ wxID_PAWLEYLOAD, wxID_PAWLEYESTIMATE, wxID_PAWLEYUPDATE,
] = [wx.NewId() for item in range(3)]

[ wxID_IMCALIBRATE,wxID_IMRECALIBRATE,wxID_IMINTEGRATE, wxID_IMCLEARCALIB,  
    wxID_IMCOPYCONTROLS, wxID_INTEGRATEALL, wxID_IMSAVECONTROLS, wxID_IMLOADCONTROLS,
] = [wx.NewId() for item in range(8)]

[ wxID_MASKCOPY, wxID_MASKSAVE, wxID_MASKLOAD,
] = [wx.NewId() for item in range(3)]

[ wxID_STRSTACOPY, wxID_STRSTAFIT, wxID_STRSTASAVE, wxID_STRSTALOAD,wxID_APPENDDZERO,
] = [wx.NewId() for item in range(5)]

[ wxID_BACKCOPY,wxID_LIMITCOPY,wxID_SAMPLECOPY, wxID_BACKFLAGCOPY, wxID_SAMPLEFLAGCOPY,
    wxID_SAMPLESAVE, wxID_SAMPLELOAD,wxID_ADDEXCLREGION,
] = [wx.NewId() for item in range(8)]

[ wxID_INSTPRMRESET,wxID_CHANGEWAVETYPE,wxID_INSTCOPY, wxID_INSTFLAGCOPY, wxID_INSTLOAD,
    wxID_INSTSAVE,
] = [wx.NewId() for item in range(6)]

[ wxID_UNDO,wxID_LSQPEAKFIT,wxID_LSQONECYCLE,wxID_RESETSIGGAM,wxID_CLEARPEAKS,wxID_AUTOSEARCH,
] = [wx.NewId() for item in range(6)]

[  wxID_INDXRELOAD, wxID_INDEXPEAKS, wxID_REFINECELL, wxID_COPYCELL, wxID_MAKENEWPHASE,
] = [wx.NewId() for item in range(5)]

[ wxID_CONSTRAINTADD,wxID_EQUIVADD,wxID_HOLDADD,wxID_FUNCTADD,
] = [wx.NewId() for item in range(4)]

[ wxID_RESTRAINTADD, wxID_RESTSELPHASE,wxID_RESTDELETE, wxID_RESRCHANGEVAL, 
    wxID_RESTCHANGEESD,wxID_AARESTRAINTADD,wxID_AARESTRAINTPLOT,
] = [wx.NewId() for item in range(7)]

[ wxID_RIGIDBODYADD,wxID_DRAWDEFINERB,wxID_RIGIDBODYIMPORT,wxID_RESIDUETORSSEQ,
    wxID_AUTOFINDRESRB,wxID_GLOBALRESREFINE,wxID_RBREMOVEALL,wxID_COPYRBPARMS,
    wxID_GLOBALTHERM,
] = [wx.NewId() for item in range(9)]

[ wxID_SAVESEQSEL,
] = [wx.NewId() for item in range(1)]

[ wxID_SELECTPHASE,
] = [wx.NewId() for item in range(1)]

[ wxID_PDFCOPYCONTROLS, wxID_PDFSAVECONTROLS, wxID_PDFLOADCONTROLS, 
    wxID_PDFCOMPUTE, wxID_PDFCOMPUTEALL, wxID_PDFADDELEMENT, wxID_PDFDELELEMENT,
] = [wx.NewId() for item in range(7)]

VERY_LIGHT_GREY = wx.Colour(235,235,235)

################################################################################
#### GSAS-II class definitions
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

      The number of keyword arguments may be increased in the future, if needs arise,
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

    :param (other): other optional keyword parameters for the
      wx.TextCtrl widget such as Size or Style may be specified.

    '''
    def __init__(self,parent,loc,key,notBlank=True,min=None,max=None,
                 OKcontrol=None,OnLeave=None,typeHint=None,
                 CIFinput=False, **kw):
        # save passed values needed outside __init__
        self.result = loc
        self.key = key
        self.OKcontrol=OKcontrol
        self.OnLeave = OnLeave
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
                self.SetValue(val)
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
                self.SetValue(val)
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
        elif val is None:
            raise Exception,("ValidatedTxtCtrl error: value of "+str(key)+
                             " element is None and typeHint not defined as int or float")
        else:
            raise Exception,("ValidatedTxtCtrl error: Unknown element ("+str(key)+
                             ") type: "+str(type(val)))
        # When the mouse is moved away or the widget loses focus
        # display the last saved value, if an expression
        self.Bind(wx.EVT_LEAVE_WINDOW, self._onLoseFocus)
        self.Bind(wx.EVT_KILL_FOCUS, self._onLoseFocus)

    def SetValue(self,val):
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
                float(val)
            except:
                if self.CIFinput and (val == '?' or val == '.'):
                    pass
                else:
                    self.invalid = True
            wx.TextCtrl.SetValue(self,str(val))
        else:
            wx.TextCtrl.SetValue(self,str(val))
            self.ShowStringValidity() # test if valid input
            return
        
        self._IndicateValidity()
        if self.OKcontrol:
            self.OKcontrol(not self.invalid)
        
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

    def _onLoseFocus(self,event):
        if self.evaluated: self.EvaluateExpression()
        if self.OnLeave: self.OnLeave(invalid=self.invalid,
                                      value=self.result[self.key],
                                      tc=self)
            
    def EvaluateExpression(self):
        '''Show the computed value when an expression is entered to the TextCtrl
        Make sure that the number fits by truncating decimal places and switching
        to scientific notation, as needed. 
        Called on loss of focus.
        '''
        if self.invalid: return # don't substitute for an invalid expression
        if not self.evaluated: return # true when an expression is evaluated
        if self.result is not None: # retrieve the stored result
            val = self.result[self.key]
            self.SetValue(G2py3.FormatValue(val))
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

################################################################################
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
class EnumSelector(wx.ComboBox):
    '''A customized :class:`wxpython.ComboBox` that selects items from a list
    of choices, but sets a dict (list) entry to the corresponding
    entry from the input list of values.

    :param wx.Panel parent: the parent to the :class:`~wxpython.ComboBox` (usually a
      frame or panel)
    :param dict dct: a dict (or list) to contain the value set
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
    :param (other): additional keyword arguments accepted by
      :class:`~wxpython.ComboBox` can be specified.
    '''
    def __init__(self,parent,dct,item,choices,values=None,**kw):
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
        self.Bind(wx.EVT_COMBOBOX, self.onSelection)
    def onSelection(self,event):
        # respond to a selection by setting the enum value in the CIF dictionary
        if self.GetValue() in self.choices: # should always be true!
            self.dct[self.item] = self.values[self.choices.index(self.GetValue())]
        else:
            self.dct[self.item] = self.values[0] # unknown

################################################################################    
def CallScrolledMultiEditor(parent,dictlst,elemlst,prelbl=[],postlbl=[],
                 title='Edit items',header='',size=(300,250),
                             CopyButton=False):
    '''Shell routine to call a ScrolledMultiEditor dialog. See
    :class:`ScrolledMultiEditor` for parameter definitions.

    :returns: True if the OK button is pressed; False if the window is closed
      with the system menu or the Cancel button.

    '''
    dlg = ScrolledMultiEditor(parent,dictlst,elemlst,prelbl,postlbl,
                              title,header,size,
                              CopyButton)
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
                 minvals=[],maxvals=[],sizevals=[]):
        if len(dictlst) != len(elemlst):
            raise Exception,"ScrolledMultiEditor error: len(dictlst) != len(elemlst) "+str(len(dictlst))+" != "+str(len(elemlst))
        wx.Dialog.__init__( # create dialog & sizer
            self,parent,wx.ID_ANY,title,
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.orig = []
        self.dictlst = dictlst
        self.elemlst = elemlst
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
        panel = wxscroll.ScrolledPanel(
            self, wx.ID_ANY,
            size=size,
            style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        cols = 3
        if CopyButton: cols += 1
        subSizer = wx.FlexGridSizer(rows=len(dictlst),cols=cols,hgap=2,vgap=2)
        self.ValidatedControlsList = [] # make list of TextCtrls
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
            if i >= len(postlbl): # label after TextCtrl, or put in a blank
                subSizer.Add((-1,-1))
            else:
                subSizer.Add(wx.StaticText(panel,wx.ID_ANY,str(postlbl[i])))
        # finish up ScrolledPanel
        panel.SetSizer(subSizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
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

    def _onClose(self,event):
        'Restore original values & close the window'
        for d,i,v in zip(self.dictlst,self.elemlst,self.orig):
            d[i] = v
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
class SymOpDialog(wx.Dialog):
    '''Class to select a symmetry operator
    '''
    def __init__(self,parent,SGData,New=True,ForceUnit=False):
        wx.Dialog.__init__(self,parent,-1,'Select symmetry operator',
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        panel = wx.Panel(self)
        self.SGData = SGData
        self.New = New
        self.Force = ForceUnit
        self.OpSelected = [0,0,0,[0,0,0],False,False]
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        if ForceUnit:
            choice = ['No','Yes']
            self.force = wx.RadioBox(panel,-1,'Force to unit cell?',choices=choice)
            self.force.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.force,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        if SGData['SGInv']:
            choice = ['No','Yes']
            self.inv = wx.RadioBox(panel,-1,'Choose inversion?',choices=choice)
            self.inv.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.inv,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        if SGData['SGLatt'] != 'P':
            LattOp = G2spc.Latt2text(SGData['SGLatt']).split(';')
            self.latt = wx.RadioBox(panel,-1,'Choose cell centering?',choices=LattOp)
            self.latt.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.latt,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        if SGData['SGLaue'] in ['-1','2/m','mmm','4/m','4/mmm']:
            Ncol = 2
        else:
            Ncol = 3
        OpList = []
        for M,T in SGData['SGOps']:
            OpList.append(G2spc.MT2text(M,T))
        self.oprs = wx.RadioBox(panel,-1,'Choose space group operator?',choices=OpList,
            majorDimension=Ncol)
        self.oprs.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
        mainSizer.Add(self.oprs,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(wx.StaticText(panel,-1,"   Choose unit cell?"),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)
        cellSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellSizer.Add((5,0),0)
        cellName = ['X','Y','Z']
        self.cell = []
        for i in range(3):
            self.cell.append(wx.SpinCtrl(panel,-1,cellName[i],size=wx.Size(50,20)))
            self.cell[-1].SetRange(-3,3)
            self.cell[-1].SetValue(0)
            self.cell[-1].Bind(wx.EVT_SPINCTRL, self.OnOpSelect)
            cellSizer.Add(self.cell[-1],0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(cellSizer,0,)
        if self.New:
            choice = ['No','Yes']
            self.new = wx.RadioBox(panel,-1,'Generate new positions?',choices=choice)
            self.new.Bind(wx.EVT_RADIOBOX, self.OnOpSelect)
            mainSizer.Add(self.new,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)

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

    def OnOpSelect(self,event):
        if self.SGData['SGInv']:
            self.OpSelected[0] = self.inv.GetSelection()
        if self.SGData['SGLatt'] != 'P':
            self.OpSelected[1] = self.latt.GetSelection()
        self.OpSelected[2] = self.oprs.GetSelection()
        for i in range(3):
            self.OpSelected[3][i] = float(self.cell[i].GetValue())
        if self.New:
            self.OpSelected[4] = self.new.GetSelection()
        if self.Force:
            self.OpSelected[5] = self.force.GetSelection()

    def GetSelection(self):
        return self.OpSelected

    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)

    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

class DisAglDialog(wx.Dialog):
    '''Distance Angle Controls dialog
    '''
    def __default__(self,data,default):
        if data:
            self.data = data
        else:
            self.data = {}
            self.data['Name'] = default['Name']
            self.data['Factors'] = [0.85,0.85]
            self.data['AtomTypes'] = default['AtomTypes']
            self.data['BondRadii'] = default['BondRadii'][:]
            self.data['AngleRadii'] = default['AngleRadii'][:]
        
    def __init__(self,parent,data,default):
        wx.Dialog.__init__(self,parent,-1,'Distance Angle Controls', 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.default = default
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.__default__(data,self.default)
        self.Draw(self.data)
                
    def Draw(self,data):
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,'Controls for phase '+data['Name']),
            0,wx.ALIGN_CENTER_VERTICAL|wx.LEFT,10)
        mainSizer.Add((10,10),1)
        
        radiiSizer = wx.FlexGridSizer(2,3,5,5)
        radiiSizer.Add(wx.StaticText(self.panel,-1,' Type'),0,wx.ALIGN_CENTER_VERTICAL)
        radiiSizer.Add(wx.StaticText(self.panel,-1,'Bond radii'),0,wx.ALIGN_CENTER_VERTICAL)
        radiiSizer.Add(wx.StaticText(self.panel,-1,'Angle radii'),0,wx.ALIGN_CENTER_VERTICAL)
        self.objList = {}
        for id,item in enumerate(self.data['AtomTypes']):
            radiiSizer.Add(wx.StaticText(self.panel,-1,' '+item),0,wx.ALIGN_CENTER_VERTICAL)
            bRadii = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['BondRadii'][id]),style=wx.TE_PROCESS_ENTER)
            self.objList[bRadii.GetId()] = ['BondRadii',id]
            bRadii.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            bRadii.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            radiiSizer.Add(bRadii,0,wx.ALIGN_CENTER_VERTICAL)
            aRadii = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['AngleRadii'][id]),style=wx.TE_PROCESS_ENTER)
            self.objList[aRadii.GetId()] = ['AngleRadii',id]
            aRadii.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            aRadii.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            radiiSizer.Add(aRadii,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(radiiSizer,0,wx.EXPAND)
        factorSizer = wx.FlexGridSizer(2,2,5,5)
        Names = ['Bond','Angle']
        for i,name in enumerate(Names):
            factorSizer.Add(wx.StaticText(self.panel,-1,name+' search factor'),0,wx.ALIGN_CENTER_VERTICAL)
            bondFact = wx.TextCtrl(self.panel,-1,value='%.3f'%(data['Factors'][i]),style=wx.TE_PROCESS_ENTER)
            self.objList[bondFact.GetId()] = ['Factors',i]
            bondFact.Bind(wx.EVT_TEXT_ENTER,self.OnRadiiVal)
            bondFact.Bind(wx.EVT_KILL_FOCUS,self.OnRadiiVal)
            factorSizer.Add(bondFact)
        mainSizer.Add(factorSizer,0,wx.EXPAND)
        
        OkBtn = wx.Button(self.panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        ResetBtn = wx.Button(self.panel,-1,'Reset')
        ResetBtn.Bind(wx.EVT_BUTTON, self.OnReset)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add(ResetBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM|wx.TOP, 10)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.Fit()
    
    def OnRadiiVal(self,event):
        Obj = event.GetEventObject()
        item = self.objList[Obj.GetId()]
        try:
            self.data[item[0]][item[1]] = float(Obj.GetValue())
        except ValueError:
            pass
        Obj.SetValue("%.3f"%(self.data[item[0]][item[1]]))          #reset in case of error
        
    def GetData(self):
        return self.data
        
    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnReset(self,event):
        data = {}
        self.__default__(data,self.default)
        self.Draw(self.data)
        
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
            lineSizer.Add(choice,0,wx.ALIGN_CENTER)
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
        
class SingleFloatDialog(wx.Dialog):
    'Dialog to obtain a single float value from user'
    def __init__(self,parent,title,prompt,value,limits=[0.,1.],format='%.5g'):
        wx.Dialog.__init__(self,parent,-1,title, 
            pos=wx.DefaultPosition,style=wx.DEFAULT_DIALOG_STYLE)
        self.panel = wx.Panel(self)         #just a dummy - gets destroyed in Draw!
        self.limits = limits
        self.value = value
        self.prompt = prompt
        self.format = format
        self.Draw()
        
    def Draw(self):
        
        def OnValItem(event):
            try:
                val = float(valItem.GetValue())
                if val < self.limits[0] or val > self.limits[1]:
                    raise ValueError
            except ValueError:
                val = self.value
            self.value = val
            valItem.SetValue(self.format%(self.value))
            
        self.panel.Destroy()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,self.prompt),0,wx.ALIGN_CENTER)
        valItem = wx.TextCtrl(self.panel,-1,value=self.format%(self.value),style=wx.TE_PROCESS_ENTER)
        mainSizer.Add(valItem,0,wx.ALIGN_CENTER)
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

    def GetValue(self):
        return self.value
        
    def OnOk(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_OK)              
        
    def OnCancel(self,event):
        parent = self.GetParent()
        parent.Raise()
        self.EndModal(wx.ID_CANCEL)

################################################################################
class SingleStringDialog(wx.Dialog):
    '''Dialog to obtain a single string value from user
    
    :param wx.Frame parent: name of parent frame
    :param str title: title string for dialog
    :param str prompt: string to tell use what they are inputting
    :param str value: default input value, if any
    '''
    def __init__(self,parent,title,prompt,value='',size=(200,-1)):
        wx.Dialog.__init__(self,parent,wx.ID_ANY,title, 
                           pos=wx.DefaultPosition,
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        self.value = value
        self.prompt = prompt
        self.CenterOnParent()
        self.panel = wx.Panel(self)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(wx.StaticText(self.panel,-1,self.prompt),0,wx.ALIGN_CENTER)
        self.valItem = wx.TextCtrl(self.panel,-1,value=str(self.value),size=size)
        mainSizer.Add(self.valItem,0,wx.ALIGN_CENTER)
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
def ItemSelector(ChoiceList, ParentFrame=None,
                 title='Select an item',
                 size=None, header='Item Selector',
                 useCancel=True,multiple=False):
        ''' Provide a wx dialog to select a single item from list of choices

        :param list ChoiceList: a list of choices where one will be selected
        :param wx.Frame ParentFrame: Name of parent frame (default None)
        :param str title: heading above list of choices (default 'Select an item')
        :param wx.Size size: Size for dialog to be created (default None -- size as needed)
        :param str header: Title to place on window frame (default 'Item Selector')
        :param bool useCancel: If True (default) both the OK and Cancel buttons are offered
        :param bool multiple: If True then multiple items can be selected (default False)

        :returns: the selection index or None or a selection list if multiple is true
        '''
        if multiple:
            if useCancel:
                dlg = wx.MultiChoiceDialog(
                    ParentFrame,title, header, ChoiceList)
            else:
                dlg = wx.MultiChoiceDialog(
                    ParentFrame,title, header, ChoiceList,
                    style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.OK|wx.CENTRE)
            pass
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
                return dlg.GetSelections()
            else:
                return dlg.GetSelection()
        else:
            return None
        dlg.Destroy()

################################################################################
class GridFractionEditor(wg.PyGridCellEditor):
    '''A grid cell editor class that allows entry of values as fractions as well
    as sine and cosine values [as s() and c()]
    '''
    def __init__(self,grid):
        wg.PyGridCellEditor.__init__(self)

    def Create(self, parent, id, evtHandler):
        self._tc = wx.TextCtrl(parent, id, "")
        self._tc.SetInsertionPoint(0)
        self.SetControl(self._tc)

        if evtHandler:
            self._tc.PushEventHandler(evtHandler)

        self._tc.Bind(wx.EVT_CHAR, self.OnChar)

    def SetSize(self, rect):
        self._tc.SetDimensions(rect.x, rect.y, rect.width+2, rect.height+2,
                               wx.SIZE_ALLOW_MINUS_ONE)

    def BeginEdit(self, row, col, grid):
        self.startValue = grid.GetTable().GetValue(row, col)
        self._tc.SetValue(str(self.startValue))
        self._tc.SetInsertionPointEnd()
        self._tc.SetFocus()
        self._tc.SetSelection(0, self._tc.GetLastPosition())

    def EndEdit(self, row, col, grid):
        changed = False

        val = self._tc.GetValue().lower()
        
        if val != self.startValue:
            changed = True
            neg = False
            if '-' in val:
                neg = True
            if '/' in val and '.' not in val:
                val += '.'
            elif 's' in val and not 'sind(' in val:
                if neg:
                    val = '-sind('+val.strip('-s')+')'
                else:
                    val = 'sind('+val.strip('s')+')'
            elif 'c' in val and not 'cosd(' in val:
                if neg:
                    val = '-cosd('+val.strip('-c')+')'
                else:
                    val = 'cosd('+val.strip('c')+')'
            try:
                val = float(eval(val))
            except (SyntaxError,NameError):
                val = self.startValue
            grid.GetTable().SetValue(row, col, val) # update the table

        self.startValue = ''
        self._tc.SetValue('')
        return changed

    def Reset(self):
        self._tc.SetValue(self.startValue)
        self._tc.SetInsertionPointEnd()

    def Clone(self):
        return GridFractionEditor(grid)

    def StartingKey(self, evt):
        self.OnChar(evt)
        if evt.GetSkipped():
            self._tc.EmulateKeyPress(evt)

    def OnChar(self, evt):
        key = evt.GetKeyCode()
        if key == 15:
            return
        if key > 255:
            evt.Skip()
            return
        char = chr(key)
        if char in '.+-/0123456789cosind()':
            self._tc.WriteText(char)
        else:
            evt.Skip()

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
        else:
            ShowHelp(helpType,self.frame)

    def OnHelpAbout(self, event):
        "Display an 'About GSAS-II' box"
        global __version__
        info = wx.AboutDialogInfo()
        info.Name = 'GSAS-II'
        ver = GSASIIpath.svnGetRev()
        if ver: 
            info.Version = 'Revision '+str(ver)+' (svn), version '+__version__
        else:
            info.Version = 'Revision '+str(GSASIIpath.GetVersionNumber())+' (.py files), version '+__version__
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
class DataFrame(wx.Frame):
    '''Create the data item window and all the entries in menus used in
    that window. For Linux and windows, the menu entries are created for the
    current data item window, but in the Mac the menu is accessed from all
    windows. This means that a different menu is posted depending on which
    data item is posted. On the Mac, all the menus contain the data tree menu
    items, but additional menus are added specific to the data item. 

    Note that while the menus are created here, 
    the binding for the menus is done later in various GSASII*GUI modules,
    where the functions to be called are defined.
    '''
    def Bind(self,*args,**kwargs):
        '''Override the Bind() function: on the Mac the binding is to
        the main window, so that menus operate with any window on top.
        For other platforms, call the default wx.Frame Bind()
        '''
        if sys.platform == "darwin": # mac
            self.G2frame.Bind(*args,**kwargs)
        else:
            wx.Frame.Bind(self,*args,**kwargs)      
        
    def PrefillDataMenu(self,menu,helpType,helpLbl=None,empty=False):
        '''Create the "standard" part of data frame menus. Note that on Linux and
        Windows nothing happens here. On Mac, this menu duplicates the
        tree menu, but adds an extra help command for the data item and a separator. 
        '''
        self.datamenu = menu
        self.helpType = helpType
        self.helpLbl = helpLbl
        if sys.platform == "darwin": # mac                         
            self.G2frame.FillMainMenu(menu) # add the data tree menu items
            if not empty:
                menu.Append(wx.Menu(title=''),title='|') # add a separator
        
    def PostfillDataMenu(self,empty=False):
        '''Create the "standard" part of data frame menus. Note that on Linux and
        Windows, this is the standard help Menu. On Mac, this menu duplicates the
        tree menu, but adds an extra help command for the data item and a separator. 
        '''
        menu = self.datamenu
        helpType = self.helpType
        helpLbl = self.helpLbl
        if sys.platform == "darwin": # mac
            if not empty:
                menu.Append(wx.Menu(title=''),title='|') # add another separator
            menu.Append(AddHelp(self.G2frame,helpType=helpType, helpLbl=helpLbl),
                        title='&Help')
        else: # other
            menu.Append(menu=MyHelp(self,helpType=helpType, helpLbl=helpLbl),
                        title='&Help')

    def _init_menus(self):
        'define all GSAS-II data frame menus'

        # for use where no menu or data frame help is provided
        self.BlankMenu = wx.MenuBar()
        
        # Controls
        self.ControlsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ControlsMenu,helpType='Controls',empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Notebook
        self.DataNotebookMenu = wx.MenuBar() 
        self.PrefillDataMenu(self.DataNotebookMenu,helpType='Notebook',empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Comments
        self.DataCommentsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DataCommentsMenu,helpType='Comments',empty=True)
        self.PostfillDataMenu(empty=True)
        
        # Constraints
        self.ConstraintMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ConstraintMenu,helpType='Constraints')
        self.ConstraintEdit = wx.Menu(title='')
        self.ConstraintMenu.Append(menu=self.ConstraintEdit, title='Edit')
        self.ConstraintEdit.Append(id=wxID_HOLDADD, kind=wx.ITEM_NORMAL,text='Add hold',
            help='Add hold on a parameter value')
        self.ConstraintEdit.Append(id=wxID_EQUIVADD, kind=wx.ITEM_NORMAL,text='Add equivalence',
            help='Add equivalence between parameter values')
        self.ConstraintEdit.Append(id=wxID_CONSTRAINTADD, kind=wx.ITEM_NORMAL,text='Add constraint',
            help='Add constraint on parameter values')
        self.ConstraintEdit.Append(id=wxID_FUNCTADD, kind=wx.ITEM_NORMAL,text='Add New Var',
            help='Add variable composed of existing parameter')
        self.PostfillDataMenu()
        
        # Rigid bodies
        self.VectorRBEdit = wx.Menu(title='')
        self.VectorRBEdit.Append(id=wxID_RIGIDBODYADD, kind=wx.ITEM_NORMAL,text='Add rigid body',
            help='Add vector rigid body')
        self.ResidueRBMenu = wx.Menu(title='')
        self.ResidueRBMenu.Append(id=wxID_RIGIDBODYIMPORT, kind=wx.ITEM_NORMAL,text='Import XYZ',
            help='Import rigid body XYZ from file')
        self.ResidueRBMenu.Append(id=wxID_RESIDUETORSSEQ, kind=wx.ITEM_NORMAL,text='Define sequence',
            help='Define torsion sequence')
        self.ResidueRBMenu.Append(id=wxID_RIGIDBODYADD, kind=wx.ITEM_NORMAL,text='Import residues',
            help='Import residue rigid bodies from macro file')
            
        self.RigidBodyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RigidBodyMenu,helpType='Rigid bodies')
        self.RigidBodyMenu.Append(menu=self.VectorRBEdit, title='Edit')        
        self.PostfillDataMenu()
            
        # Restraints
        self.RestraintEdit = wx.Menu(title='')
        self.RestraintEdit.Append(id=wxID_RESTSELPHASE, kind=wx.ITEM_NORMAL,text='Select phase',
            help='Select phase')
        self.RestraintEdit.Append(id=wxID_RESTRAINTADD, kind=wx.ITEM_NORMAL,text='Add restraints',
            help='Add restraints')
        self.RestraintEdit.Enable(wxID_RESTRAINTADD,True)    #gets disenabled if macromolecule phase
        self.RestraintEdit.Append(id=wxID_AARESTRAINTADD, kind=wx.ITEM_NORMAL,text='Add residue restraints',
            help='Add residue based restraints for macromolecules from macro file')
        self.RestraintEdit.Enable(wxID_AARESTRAINTADD,False)    #gets enabled if macromolecule phase
        self.RestraintEdit.Append(id=wxID_AARESTRAINTPLOT, kind=wx.ITEM_NORMAL,text='Plot residue restraints',
            help='Plot selected residue based restraints for macromolecules from macro file')
        self.RestraintEdit.Enable(wxID_AARESTRAINTPLOT,False)    #gets enabled if macromolecule phase
        self.RestraintEdit.Append(id=wxID_RESRCHANGEVAL, kind=wx.ITEM_NORMAL,text='Change value',
            help='Change observed value')
        self.RestraintEdit.Append(id=wxID_RESTCHANGEESD, kind=wx.ITEM_NORMAL,text='Change esd',
            help='Change esd in observed value')
        self.RestraintEdit.Append(id=wxID_RESTDELETE, kind=wx.ITEM_NORMAL,text='Delete restraints',
            help='Delete selected restraints')

        self.RestraintMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RestraintMenu,helpType='Restraints')
        self.RestraintMenu.Append(menu=self.RestraintEdit, title='Edit')
        self.PostfillDataMenu()
            
        # Sequential results
        self.SequentialMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SequentialMenu,helpType='Sequential',helpLbl='Sequential Refinement')
        self.SequentialFile = wx.Menu(title='')
        self.SequentialMenu.Append(menu=self.SequentialFile, title='File')
        self.SequentialFile.Append(id=wxID_SAVESEQSEL, kind=wx.ITEM_NORMAL,text='Save...',
            help='Save selected sequential refinement results')
        self.PostfillDataMenu()
            
        # PDR
        self.ErrorMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ErrorMenu,helpType='PWD Analysis',helpLbl='Powder Fit Error Analysis')
        self.ErrorAnal = wx.Menu(title='')
        self.ErrorMenu.Append(menu=self.ErrorAnal,title='Analysis')
        self.ErrorAnal.Append(id=wxID_PWDANALYSIS,kind=wx.ITEM_NORMAL,text='Analyze',
            help='Error analysis on powder pattern')
        self.PostfillDataMenu()
            
        # PDR / Limits
        self.LimitMenu = wx.MenuBar()
        self.PrefillDataMenu(self.LimitMenu,helpType='Limits')
        self.LimitEdit = wx.Menu(title='')
        self.LimitMenu.Append(menu=self.LimitEdit, title='Edit')
        self.LimitEdit.Append(id=wxID_LIMITCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy limits to other histograms')
        self.LimitEdit.Append(id=wxID_ADDEXCLREGION, kind=wx.ITEM_NORMAL,text='Add exclude',
            help='Add excluded region - select a point on plot; drag to adjust')            
        self.PostfillDataMenu()
            
        # PDR / Background
        self.BackMenu = wx.MenuBar()
        self.PrefillDataMenu(self.BackMenu,helpType='Background')
        self.BackEdit = wx.Menu(title='')
        self.BackMenu.Append(menu=self.BackEdit, title='File')
        self.BackEdit.Append(id=wxID_BACKCOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy background parameters to other histograms')
        self.BackEdit.Append(id=wxID_BACKFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy background refinement flags to other histograms')
        self.PostfillDataMenu()
            
        # PDR / Instrument Parameters
        self.InstMenu = wx.MenuBar()
        self.PrefillDataMenu(self.InstMenu,helpType='Instrument Parameters')
        self.InstEdit = wx.Menu(title='')
        self.InstMenu.Append(menu=self.InstEdit, title='Operations')
        self.InstEdit.Append(help='Reset instrument profile parameters to default', 
            id=wxID_INSTLOAD, kind=wx.ITEM_NORMAL,text='Load profile...')
        self.InstEdit.Append(help='Load instrument profile parameters from file', 
            id=wxID_INSTSAVE, kind=wx.ITEM_NORMAL,text='Save profile...')
        self.InstEdit.Append(help='Save instrument profile parameters to file', 
            id=wxID_INSTPRMRESET, kind=wx.ITEM_NORMAL,text='Reset profile')
        self.InstEdit.Append(help='Copy instrument profile parameters to other histograms', 
            id=wxID_INSTCOPY, kind=wx.ITEM_NORMAL,text='Copy')
        self.InstEdit.Append(id=wxID_INSTFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy instrument parameter refinement flags to other histograms')
#        self.InstEdit.Append(help='Change radiation type (Ka12 - synch)', 
#            id=wxID_CHANGEWAVETYPE, kind=wx.ITEM_NORMAL,text='Change radiation')
        self.PostfillDataMenu()
        
        # PDR / Sample Parameters
        self.SampleMenu = wx.MenuBar()
        self.PrefillDataMenu(self.SampleMenu,helpType='Sample Parameters')
        self.SampleEdit = wx.Menu(title='')
        self.SampleMenu.Append(menu=self.SampleEdit, title='File')
        self.SampleEdit.Append(id=wxID_SAMPLELOAD, kind=wx.ITEM_NORMAL,text='Load',
            help='Load sample parameters from file')
        self.SampleEdit.Append(id=wxID_SAMPLESAVE, kind=wx.ITEM_NORMAL,text='Save',
            help='Save sample parameters to file')
        self.SampleEdit.Append(id=wxID_SAMPLECOPY, kind=wx.ITEM_NORMAL,text='Copy',
            help='Copy refinable sample parameters to other histograms')
        self.SampleEdit.Append(id=wxID_SAMPLEFLAGCOPY, kind=wx.ITEM_NORMAL,text='Copy flags',
            help='Copy sample parameter refinement flags to other histograms')
        self.PostfillDataMenu()

        # PDR / Peak List
        self.PeakMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PeakMenu,helpType='Peak List')
        self.PeakEdit = wx.Menu(title='')
        self.PeakMenu.Append(menu=self.PeakEdit, title='Peak Fitting')
        self.AutoSearch = self.PeakEdit.Append(help='Automatic peak search', 
            id=wxID_AUTOSEARCH, kind=wx.ITEM_NORMAL,text='Auto search')
        self.UnDo = self.PeakEdit.Append(help='Undo last least squares refinement', 
            id=wxID_UNDO, kind=wx.ITEM_NORMAL,text='UnDo')
        self.PeakFit = self.PeakEdit.Append(id=wxID_LSQPEAKFIT, kind=wx.ITEM_NORMAL,text='LSQ PeakFit', 
            help='Peak fitting via least-squares' )
        self.PFOneCycle = self.PeakEdit.Append(id=wxID_LSQONECYCLE, kind=wx.ITEM_NORMAL,text='LSQ one cycle', 
            help='One cycle of Peak fitting via least-squares' )
        self.PeakEdit.Append(id=wxID_RESETSIGGAM, kind=wx.ITEM_NORMAL, 
            text='Reset sig and gam',help='Reset sigma and gamma to global fit' )
        self.PeakEdit.Append(id=wxID_CLEARPEAKS, kind=wx.ITEM_NORMAL,text='Clear peaks', 
            help='Clear the peak list' )
        self.PostfillDataMenu()
        self.UnDo.Enable(False)
        self.PeakFit.Enable(False)
        self.PFOneCycle.Enable(False)
        
        # PDR / Index Peak List
        self.IndPeaksMenu = wx.MenuBar()
        self.PrefillDataMenu(self.IndPeaksMenu,helpType='Index Peak List')
        self.IndPeaksEdit = wx.Menu(title='')
        self.IndPeaksMenu.Append(menu=self.IndPeaksEdit,title='Operations')
        self.IndPeaksEdit.Append(help='Load/Reload index peaks from peak list',id=wxID_INDXRELOAD, 
            kind=wx.ITEM_NORMAL,text='Load/Reload')
        self.PostfillDataMenu()
        
        # PDR / Unit Cells List
        self.IndexMenu = wx.MenuBar()
        self.PrefillDataMenu(self.IndexMenu,helpType='Unit Cells List')
        self.IndexEdit = wx.Menu(title='')
        self.IndexMenu.Append(menu=self.IndexEdit, title='Cell Index/Refine')
        self.IndexPeaks = self.IndexEdit.Append(help='', id=wxID_INDEXPEAKS, kind=wx.ITEM_NORMAL,
            text='Index Cell')
        self.CopyCell = self.IndexEdit.Append( id=wxID_COPYCELL, kind=wx.ITEM_NORMAL,text='Copy Cell', 
            help='Copy selected unit cell from indexing to cell refinement fields')
        self.RefineCell = self.IndexEdit.Append( id=wxID_REFINECELL, kind=wx.ITEM_NORMAL, 
            text='Refine Cell',help='Refine unit cell parameters from indexed peaks')
        self.MakeNewPhase = self.IndexEdit.Append( id=wxID_MAKENEWPHASE, kind=wx.ITEM_NORMAL,
            text='Make new phase',help='Make new phase from selected unit cell')
        self.PostfillDataMenu()
        self.IndexPeaks.Enable(False)
        self.CopyCell.Enable(False)
        self.RefineCell.Enable(False)
        self.MakeNewPhase.Enable(False)
        
        # PDR / Reflection Lists
        self.ReflMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ReflMenu,helpType='Reflection List')
        self.ReflEdit = wx.Menu(title='')
        self.ReflMenu.Append(menu=self.ReflEdit, title='Reflection List')
        self.SelectPhase = self.ReflEdit.Append(help='Select phase for reflection list',id=wxID_SELECTPHASE, 
            kind=wx.ITEM_NORMAL,text='Select phase')
        self.PostfillDataMenu()
        
        # IMG / Image Controls
        self.ImageMenu = wx.MenuBar()
        self.PrefillDataMenu(self.ImageMenu,helpType='Image Controls')
        self.ImageEdit = wx.Menu(title='')
        self.ImageMenu.Append(menu=self.ImageEdit, title='Operations')
        self.ImageEdit.Append(help='Calibrate detector by fitting to calibrant lines', 
            id=wxID_IMCALIBRATE, kind=wx.ITEM_NORMAL,text='Calibrate')
        self.ImageEdit.Append(help='Recalibrate detector by fitting to calibrant lines', 
            id=wxID_IMRECALIBRATE, kind=wx.ITEM_NORMAL,text='Recalibrate')
        self.ImageEdit.Append(help='Clear calibration data points and rings',id=wxID_IMCLEARCALIB, 
            kind=wx.ITEM_NORMAL,text='Clear calibration')
        self.ImageEdit.Append(help='Integrate selected image',id=wxID_IMINTEGRATE, 
            kind=wx.ITEM_NORMAL,text='Integrate')
        self.ImageEdit.Append(help='Integrate all images selected from list',id=wxID_INTEGRATEALL,
            kind=wx.ITEM_NORMAL,text='Integrate all')
        self.ImageEdit.Append(help='Copy image controls to other images', 
            id=wxID_IMCOPYCONTROLS, kind=wx.ITEM_NORMAL,text='Copy Controls')
        self.ImageEdit.Append(help='Save image controls to file', 
            id=wxID_IMSAVECONTROLS, kind=wx.ITEM_NORMAL,text='Save Controls')
        self.ImageEdit.Append(help='Load image controls from file', 
            id=wxID_IMLOADCONTROLS, kind=wx.ITEM_NORMAL,text='Load Controls')
        self.PostfillDataMenu()
            
        # IMG / Masks
        self.MaskMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MaskMenu,helpType='Image Masks')
        self.MaskEdit = wx.Menu(title='')
        self.MaskMenu.Append(menu=self.MaskEdit, title='Operations')
        self.MaskEdit.Append(help='Copy mask to other images', 
            id=wxID_MASKCOPY, kind=wx.ITEM_NORMAL,text='Copy mask')
        self.MaskEdit.Append(help='Save mask to file', 
            id=wxID_MASKSAVE, kind=wx.ITEM_NORMAL,text='Save mask')
        self.MaskEdit.Append(help='Load mask from file', 
            id=wxID_MASKLOAD, kind=wx.ITEM_NORMAL,text='Load mask')
        self.PostfillDataMenu()
            
        # IMG / Stress/Strain
        self.StrStaMenu = wx.MenuBar()
        self.PrefillDataMenu(self.StrStaMenu,helpType='Stress/Strain')
        self.StrStaEdit = wx.Menu(title='')
        self.StrStaMenu.Append(menu=self.StrStaEdit, title='Operations')
        self.StrStaEdit.Append(help='Append d-zero for one ring', 
            id=wxID_APPENDDZERO, kind=wx.ITEM_NORMAL,text='Append d-zero')
        self.StrStaEdit.Append(help='Fit stress/strain data', 
            id=wxID_STRSTAFIT, kind=wx.ITEM_NORMAL,text='Fit stress/strain')
        self.StrStaEdit.Append(help='Copy stress/strain data to other images', 
            id=wxID_STRSTACOPY, kind=wx.ITEM_NORMAL,text='Copy stress/strain')
        self.StrStaEdit.Append(help='Save stress/strain data to file', 
            id=wxID_STRSTASAVE, kind=wx.ITEM_NORMAL,text='Save stress/strain')
        self.StrStaEdit.Append(help='Load stress/strain data from file', 
            id=wxID_STRSTALOAD, kind=wx.ITEM_NORMAL,text='Load stress/strain')
        self.PostfillDataMenu()
            
        # PDF / PDF Controls
        self.PDFMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PDFMenu,helpType='PDF Controls')
        self.PDFEdit = wx.Menu(title='')
        self.PDFMenu.Append(menu=self.PDFEdit, title='PDF Controls')
        self.PDFEdit.Append(help='Add element to sample composition',id=wxID_PDFADDELEMENT, kind=wx.ITEM_NORMAL,
            text='Add element')
        self.PDFEdit.Append(help='Delete element from sample composition',id=wxID_PDFDELELEMENT, kind=wx.ITEM_NORMAL,
            text='Delete element')
        self.PDFEdit.Append(help='Copy PDF controls', id=wxID_PDFCOPYCONTROLS, kind=wx.ITEM_NORMAL,
            text='Copy controls')
        #        self.PDFEdit.Append(help='Load PDF controls from file',id=wxID_PDFLOADCONTROLS, kind=wx.ITEM_NORMAL,
        #            text='Load Controls')
        #        self.PDFEdit.Append(help='Save PDF controls to file', id=wxID_PDFSAVECONTROLS, kind=wx.ITEM_NORMAL,
        #            text='Save controls')
        self.PDFEdit.Append(help='Compute PDF', id=wxID_PDFCOMPUTE, kind=wx.ITEM_NORMAL,
            text='Compute PDF')
        self.PDFEdit.Append(help='Compute all PDFs', id=wxID_PDFCOMPUTEALL, kind=wx.ITEM_NORMAL,
            text='Compute all PDFs')
        self.PostfillDataMenu()
        
        # Phase / General tab
        self.DataGeneral = wx.MenuBar()
        self.PrefillDataMenu(self.DataGeneral,helpType='General', helpLbl='Phase/General')
        self.DataGeneral.Append(menu=wx.Menu(title=''),title='Select tab')
        self.GeneralCalc = wx.Menu(title='')
        self.DataGeneral.Append(menu=self.GeneralCalc,title='Compute')
        self.GeneralCalc.Append(help='Compute Fourier map',id=wxID_FOURCALC, kind=wx.ITEM_NORMAL,
            text='Fourier map')
        self.GeneralCalc.Append(help='Search Fourier map',id=wxID_FOURSEARCH, kind=wx.ITEM_NORMAL,
            text='Search map')
        self.GeneralCalc.Append(help='Run charge flipping',id=wxID_CHARGEFLIP, kind=wx.ITEM_NORMAL,
            text='Charge flipping')
        self.GeneralCalc.Append(help='Clear map',id=wxID_FOURCLEAR, kind=wx.ITEM_NORMAL,
            text='Clear map')
        self.GeneralCalc.Append(help='Run Monte Carlo - Simulated Annealing',id=wxID_SINGLEMCSA, kind=wx.ITEM_NORMAL,
            text='MC/SA')
        self.GeneralCalc.Append(help='Run Monte Carlo - Simulated Annealing on multiprocessors',id=wxID_MULTIMCSA, kind=wx.ITEM_NORMAL,
            text='Multi MC/SA')            #currently not useful
        self.PostfillDataMenu()
        
        # Phase / Data tab
        self.DataMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DataMenu,helpType='Data', helpLbl='Phase/Data')
        self.DataMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.DataEdit = wx.Menu(title='')
        self.DataMenu.Append(menu=self.DataEdit, title='Edit')
        self.DataEdit.Append(id=wxID_PWDRADD, kind=wx.ITEM_NORMAL,text='Add powder histograms',
            help='Select new powder histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_HKLFADD, kind=wx.ITEM_NORMAL,text='Add single crystal histograms',
            help='Select new single crystal histograms to be used for this phase')
        self.DataEdit.Append(id=wxID_DATADELETE, kind=wx.ITEM_NORMAL,text='Delete histograms',
            help='Delete histograms from use for this phase')
        self.PostfillDataMenu()
            
        # Phase / Atoms tab
        self.AtomsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.AtomsMenu,helpType='Atoms')
        self.AtomsMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.AtomEdit = wx.Menu(title='')
        self.AtomCompute = wx.Menu(title='')
        self.AtomsMenu.Append(menu=self.AtomEdit, title='Edit')
        self.AtomsMenu.Append(menu=self.AtomCompute, title='Compute')
        self.AtomEdit.Append(id=wxID_ATOMSEDITADD, kind=wx.ITEM_NORMAL,text='Append atom',
            help='Appended as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSVIEWADD, kind=wx.ITEM_NORMAL,text='Append view point',
            help='Appended as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMSEDITINSERT, kind=wx.ITEM_NORMAL,text='Insert atom',
            help='Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMVIEWINSERT, kind=wx.ITEM_NORMAL,text='Insert view point',
            help='Select atom row to insert before; inserted as an H atom')
        self.AtomEdit.Append(id=wxID_ATOMMOVE, kind=wx.ITEM_NORMAL,text='Move atom to view point',
            help='Select single atom to move')
        self.AtomEdit.Append(id=wxID_ATOMSEDITDELETE, kind=wx.ITEM_NORMAL,text='Delete atom',
            help='Select atoms to delete first')
        self.AtomEdit.Append(id=wxID_ATOMSREFINE, kind=wx.ITEM_NORMAL,text='Set atom refinement flags',
            help='Select atoms to refine first')
        self.AtomEdit.Append(id=wxID_ATOMSMODIFY, kind=wx.ITEM_NORMAL,text='Modify atom parameters',
            help='Select atoms to modify first')
        self.AtomEdit.Append(id=wxID_ATOMSTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform atoms',
            help='Select atoms to transform first')
        self.AtomEdit.Append(id=wxID_RELOADDRAWATOMS, kind=wx.ITEM_NORMAL,text='Reload draw atoms',
            help='Reload atom drawing list')
        submenu = wx.Menu()
        self.AtomEdit.AppendMenu(wx.ID_ANY, 'Reimport atoms', submenu, 
            help='Reimport atoms from file; sequence must match')
        # setup a cascade menu for the formats that have been defined
        self.ReImportMenuId = {}  # points to readers for each menu entry
        for reader in self.G2frame.ImportPhaseReaderlist:
            item = submenu.Append(
                wx.ID_ANY,help=reader.longFormatName,
                kind=wx.ITEM_NORMAL,text='reimport coordinates from '+reader.formatName+' file')
            self.ReImportMenuId[item.GetId()] = reader
        item = submenu.Append(
            wx.ID_ANY,
            help='Reimport coordinates, try to determine format from file',
            kind=wx.ITEM_NORMAL,
            text='guess format from file')
        self.ReImportMenuId[item.GetId()] = None # try all readers

        self.AtomCompute.Append(id=wxID_ATOMSDISAGL, kind=wx.ITEM_NORMAL,text='Show Distances && Angles',
            help='Compute distances & angles for selected atoms')
        self.AtomCompute.Append(id=wxID_ATOMSPDISAGL, kind=wx.ITEM_NORMAL,text='Save Distances && Angles',
            help='Compute distances & angles for selected atoms')
        self.PostfillDataMenu()
                 
        # Phase / Draw Options tab
        self.DataDrawOptions = wx.MenuBar()
        self.PrefillDataMenu(self.DataDrawOptions,helpType='Draw Options', helpLbl='Phase/Draw Options')
        self.DataDrawOptions.Append(menu=wx.Menu(title=''),title='Select tab')
        self.PostfillDataMenu()
        
        # Phase / Draw Atoms tab 
        self.DrawAtomsMenu = wx.MenuBar()
        self.PrefillDataMenu(self.DrawAtomsMenu,helpType='Draw Atoms')
        self.DrawAtomsMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.DrawAtomEdit = wx.Menu(title='')
        self.DrawAtomCompute = wx.Menu(title='')
        self.DrawAtomRestraint = wx.Menu(title='')
        self.DrawAtomRigidBody = wx.Menu(title='')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomEdit, title='Edit')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomCompute,title='Compute')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomRestraint, title='Restraints')
        self.DrawAtomsMenu.Append(menu=self.DrawAtomRigidBody, title='Rigid body')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMSTYLE, kind=wx.ITEM_NORMAL,text='Atom style',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMLABEL, kind=wx.ITEM_NORMAL,text='Atom label',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMCOLOR, kind=wx.ITEM_NORMAL,text='Atom color',
            help='Select atoms first')
        self.DrawAtomEdit.Append(id=wxID_DRAWATOMRESETCOLOR, kind=wx.ITEM_NORMAL,text='Reset atom colors',
            help='Resets all atom colors to defaults')
        self.DrawAtomEdit.Append(id=wxID_DRAWVIEWPOINT, kind=wx.ITEM_NORMAL,text='View point',
            help='View point is 1st atom selected')
        self.DrawAtomEdit.Append(id=wxID_DRAWADDEQUIV, kind=wx.ITEM_NORMAL,text='Add atoms',
            help='Add symmetry & cell equivalents to drawing set from selected atoms')
        self.DrawAtomEdit.Append(id=wxID_DRAWTRANSFORM, kind=wx.ITEM_NORMAL,text='Transform atoms',
            help='Transform selected atoms by symmetry & cell translations')
        self.DrawAtomEdit.Append(id=wxID_DRAWFILLCOORD, kind=wx.ITEM_NORMAL,text='Fill CN-sphere',
            help='Fill coordination sphere for selected atoms')            
        self.DrawAtomEdit.Append(id=wxID_DRAWFILLCELL, kind=wx.ITEM_NORMAL,text='Fill unit cell',
            help='Fill unit cell with selected atoms')
        self.DrawAtomEdit.Append(id=wxID_DRAWDELETE, kind=wx.ITEM_NORMAL,text='Delete atoms',
            help='Delete atoms from drawing set')
        self.DrawAtomCompute.Append(id=wxID_DRAWDISTVP, kind=wx.ITEM_NORMAL,text='View pt. dist.',
            help='Compute distance of selected atoms from view point')   
        self.DrawAtomCompute.Append(id=wxID_DRAWDISAGLTOR, kind=wx.ITEM_NORMAL,text='Dist. Ang. Tors.',
            help='Compute distance, angle or torsion for 2-4 selected atoms')   
        self.DrawAtomCompute.Append(id=wxID_DRAWPLANE, kind=wx.ITEM_NORMAL,text='Best plane',
            help='Compute best plane for 4+ selected atoms')   
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRBOND, kind=wx.ITEM_NORMAL,text='Add bond restraint',
            help='Add bond restraint for selected atoms (2)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRANGLE, kind=wx.ITEM_NORMAL,text='Add angle restraint',
            help='Add angle restraint for selected atoms (3: one end 1st)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRPLANE, kind=wx.ITEM_NORMAL,text='Add plane restraint',
            help='Add plane restraint for selected atoms (4+)')
        self.DrawAtomRestraint.Append(id=wxID_DRAWRESTRCHIRAL, kind=wx.ITEM_NORMAL,text='Add chiral restraint',
            help='Add chiral restraint for selected atoms (4: center atom 1st)')
        self.DrawAtomRigidBody.Append(id=wxID_DRAWDEFINERB, kind=wx.ITEM_NORMAL,text='Define rigid body',
            help='Define rigid body with selected atoms')
        self.PostfillDataMenu()

        # Phase / MCSA tab
        self.MCSAMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MCSAMenu,helpType='MC/SA')
        self.MCSAMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.MCSAEdit = wx.Menu(title='')
        self.MCSAMenu.Append(menu=self.MCSAEdit, title='MC/SA')
        self.MCSAEdit.Append(id=wxID_ADDMCSAATOM, kind=wx.ITEM_NORMAL,text='Add atom', 
            help='Add single atom to MC/SA model')
        self.MCSAEdit.Append(id=wxID_ADDMCSARB, kind=wx.ITEM_NORMAL,text='Add rigid body', 
            help='Add rigid body to MC/SA model' )
        self.MCSAEdit.Append(id=wxID_CLEARMCSARB, kind=wx.ITEM_NORMAL,text='Clear rigid bodies', 
            help='Clear all atoms & rigid bodies from MC/SA model' )
        self.MCSAEdit.Append(id=wxID_MOVEMCSA, kind=wx.ITEM_NORMAL,text='Move MC/SA solution', 
            help='Move MC/SA solution to atom list' )
        self.MCSAEdit.Append(id=wxID_MCSACLEARRESULTS, kind=wx.ITEM_NORMAL,text='Clear results', 
            help='Clear table of MC/SA results' )
        self.PostfillDataMenu()
            
        # Phase / Texture tab
        self.TextureMenu = wx.MenuBar()
        self.PrefillDataMenu(self.TextureMenu,helpType='Texture')
        self.TextureMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.TextureEdit = wx.Menu(title='')
        self.TextureMenu.Append(menu=self.TextureEdit, title='Texture')
        self.TextureEdit.Append(id=wxID_REFINETEXTURE, kind=wx.ITEM_NORMAL,text='Refine texture', 
            help='Refine the texture coefficients from sequential Pawley results')
        self.TextureEdit.Append(id=wxID_CLEARTEXTURE, kind=wx.ITEM_NORMAL,text='Clear texture', 
            help='Clear the texture coefficients' )
        self.PostfillDataMenu()
            
        # Phase / Pawley tab
        self.PawleyMenu = wx.MenuBar()
        self.PrefillDataMenu(self.PawleyMenu,helpType='Pawley')
        self.PawleyMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.PawleyEdit = wx.Menu(title='')
        self.PawleyMenu.Append(menu=self.PawleyEdit,title='Operations')
        self.PawleyEdit.Append(id=wxID_PAWLEYLOAD, kind=wx.ITEM_NORMAL,text='Pawley create',
            help='Initialize Pawley reflection list')
        self.PawleyEdit.Append(id=wxID_PAWLEYESTIMATE, kind=wx.ITEM_NORMAL,text='Pawley estimate',
            help='Estimate initial Pawley intensities')
        self.PawleyEdit.Append(id=wxID_PAWLEYUPDATE, kind=wx.ITEM_NORMAL,text='Pawley update',
            help='Update negative Pawley intensities with -0.5*Fobs and turn off refinemnt')
        self.PostfillDataMenu()
            
        # Phase / Map peaks tab
        self.MapPeaksMenu = wx.MenuBar()
        self.PrefillDataMenu(self.MapPeaksMenu,helpType='Map peaks')
        self.MapPeaksMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.MapPeaksEdit = wx.Menu(title='')
        self.MapPeaksMenu.Append(menu=self.MapPeaksEdit, title='Map peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSMOVE, kind=wx.ITEM_NORMAL,text='Move peaks', 
            help='Move selected peaks to atom list')
        self.MapPeaksEdit.Append(id=wxID_PEAKSVIEWPT, kind=wx.ITEM_NORMAL,text='View point',
            help='View point is 1st peak selected')
        self.MapPeaksEdit.Append(id=wxID_PEAKSDISTVP, kind=wx.ITEM_NORMAL,text='View pt. dist.',
            help='Compute distance of selected peaks from view point')   
        self.MapPeaksEdit.Append(id=wxID_SHOWBONDS, kind=wx.ITEM_NORMAL,text='Hide bonds',
            help='Hide or show bonds between peak positions')   
        self.MapPeaksEdit.Append(id=wxID_PEAKSDA, kind=wx.ITEM_NORMAL,text='Calc dist/ang', 
            help='Calculate distance or angle for selection')
        self.MapPeaksEdit.Append(id=wxID_FINDEQVPEAKS, kind=wx.ITEM_NORMAL,text='Equivalent peaks', 
            help='Find equivalent peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSUNIQUE, kind=wx.ITEM_NORMAL,text='Unique peaks', 
            help='Select unique set')
        self.MapPeaksEdit.Append(id=wxID_PEAKSDELETE, kind=wx.ITEM_NORMAL,text='Delete peaks', 
            help='Delete selected peaks')
        self.MapPeaksEdit.Append(id=wxID_PEAKSCLEAR, kind=wx.ITEM_NORMAL,text='Clear peaks', 
            help='Clear the map peak list')
        self.PostfillDataMenu()

        # Phase / Rigid bodies tab
        self.RigidBodiesMenu = wx.MenuBar()
        self.PrefillDataMenu(self.RigidBodiesMenu,helpType='Rigid bodies')
        self.RigidBodiesMenu.Append(menu=wx.Menu(title=''),title='Select tab')
        self.RigidBodiesEdit = wx.Menu(title='')
        self.RigidBodiesMenu.Append(menu=self.RigidBodiesEdit, title='Edit')
        self.RigidBodiesEdit.Append(id=wxID_ASSIGNATMS2RB, kind=wx.ITEM_NORMAL,text='Assign atoms to rigid body',
            help='Select & position rigid body in structure of existing atoms')
        self.RigidBodiesEdit.Append(id=wxID_AUTOFINDRESRB, kind=wx.ITEM_NORMAL,text='Auto find residues',
            help='Auto find of residue RBs in macromolecule')
        self.RigidBodiesEdit.Append(id=wxID_COPYRBPARMS, kind=wx.ITEM_NORMAL,text='Copy rigid body parms',
            help='Copy rigid body location & TLS parameters')
        self.RigidBodiesEdit.Append(id=wxID_GLOBALTHERM, kind=wx.ITEM_NORMAL,text='Global thermal motion',
            help='Global setting of residue thermal motion models')
        self.RigidBodiesEdit.Append(id=wxID_GLOBALRESREFINE, kind=wx.ITEM_NORMAL,text='Global residue refine',
            help='Global setting of residue RB refinement flags')
        self.RigidBodiesEdit.Append(id=wxID_RBREMOVEALL, kind=wx.ITEM_NORMAL,text='Remove all rigid bodies',
            help='Remove all rigid body assignment for atoms')
        self.PostfillDataMenu()
    # end of GSAS-II menu definitions
        
    def _init_ctrls(self, parent,name=None,size=None,pos=None):
        wx.Frame.__init__(self,parent=parent,
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX | wx.FRAME_FLOAT_ON_PARENT ,
            size=size,pos=pos,title='GSAS-II data display')
        self._init_menus()
        if name:
            self.SetLabel(name)
        self.Show()
        
    def __init__(self,parent,frame,data=None,name=None, size=None,pos=None):
        self.G2frame = frame
        self._init_ctrls(parent,name,size,pos)
        self.data = data
        clientSize = wx.ClientDisplayRect()
        Size = self.GetSize()
        xPos = clientSize[2]-Size[0]
        self.SetPosition(wx.Point(xPos,clientSize[1]+250))
        self.AtomGrid = []
        self.selectedRow = 0
        
    def setSizePosLeft(self,Width):
        clientSize = wx.ClientDisplayRect()
        Width[1] = min(Width[1],clientSize[2]-300)
        Width[0] = max(Width[0],300)
        self.SetSize(Width)
#        self.SetPosition(wx.Point(clientSize[2]-Width[0],clientSize[1]+250))
        
    def Clear(self):
        self.ClearBackground()
        self.DestroyChildren()
                   
################################################################################
#####  GSNotebook
################################################################################           
       
class GSNoteBook(wx.aui.AuiNotebook):
    '''Notebook used in various locations; implemented with wx.aui extension
    '''
    def __init__(self, parent, name='',size = None):
        wx.aui.AuiNotebook.__init__(self, parent, -1,
                                    style=wx.aui.AUI_NB_TOP |
                                    wx.aui.AUI_NB_SCROLL_BUTTONS)
        if size: self.SetSize(size)
                                                      
    def Clear(self):        
        GSNoteBook.DeleteAllPages(self)
        
    def FindPage(self,name):
        numPage = self.GetPageCount()
        for page in range(numPage):
            if self.GetPageText(page) == name:
                return page

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
            
################################################################################
#####  GSGrid
################################################################################           
       
class GSGrid(wg.Grid):
    '''Basic wx.Grid implementation
    '''
    def __init__(self, parent, name=''):
        wg.Grid.__init__(self,parent,-1,name=name)                    
        #self.SetSize(parent.GetClientSize())
        # above removed to speed drawing of initial grid
        # does not appear to be needed
            
    def Clear(self):
        wg.Grid.ClearGrid(self)
        
    def SetCellStyle(self,r,c,color="white",readonly=True):
        self.SetCellBackgroundColour(r,c,color)
        self.SetReadOnly(r,c,isReadOnly=readonly)
        
    def GetSelection(self):
        #this is to satisfy structure drawing stuff in G2plt when focus changes
        return None
                                                
################################################################################
#####  Table
################################################################################           
       
class Table(wg.PyGridTableBase):
    '''Basic data table for use with GSgrid
    '''
    def __init__(self, data=[], rowLabels=None, colLabels=None, types = None):
        wg.PyGridTableBase.__init__(self)
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
            if irow <> pos:
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
            return self.dataTypes[col]
        except TypeError:
            return None

    def GetValue(self, row, col):
        try:
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
            except IndexError:
                print row,col,value
                # add a new row
                if row > self.GetNumberRows():
                    self.data.append([''] * self.GetNumberCols())
                elif col > self.GetNumberCols():
                    for row in range(self.GetNumberRows):
                        self.data[row].append('')
                print self.data
                self.data[row][col] = value
        innerSetValue(row, col, value)
        
################################################################################
#### Help
################################################################################

def ShowHelp(helpType,frame):
    '''Called to bring up a web page for documentation.'''
    global htmlFirstUse
    # look up a definition for help info from dict
    helplink = helpLocDict.get(helpType)
    if helplink is None:
        # no defined link to use, create a default based on key
        helplink = 'gsasII.html#'+helpType.replace(' ','_')
    helplink = os.path.join(path2GSAS2,'help',helplink)
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

################################################################################
#####  Notebook
################################################################################           
       
def UpdateNotebook(G2frame,data):
    '''Called when the data tree notebook entry is selected. Allows for
    editing of the text in that tree entry
    '''
    def OnNoteBook(event):
        data = G2frame.dataDisplay.GetValue()
        G2frame.PatternTree.SetItemPyData(GetPatternTreeItemId(G2frame,G2frame.root,'Notebook'),data)
                    
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.dataFrame.SetLabel('Notebook')
    G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
        style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER | wx.TE_DONTWRAP)
    G2frame.dataDisplay.Bind(wx.EVT_TEXT_ENTER,OnNoteBook)
    G2frame.dataDisplay.Bind(wx.EVT_KILL_FOCUS,OnNoteBook)
    for line in data:
        G2frame.dataDisplay.AppendText(line+"\n")
    G2frame.dataDisplay.AppendText('Notebook entry @ '+time.ctime()+"\n")
    G2frame.dataFrame.setSizePosLeft([400,250])
            
################################################################################
#####  Controls
################################################################################           
       
def UpdateControls(G2frame,data):
    '''Edit overall GSAS-II controls in main Controls data tree entry
    '''
    #patch
    if 'deriv type' not in data:
        data = {}
        data['deriv type'] = 'analytic Hessian'
        data['min dM/M'] = 0.0001
        data['shift factor'] = 1.
        data['max cyc'] = 3        
        data['F**2'] = True
        data['minF/sig'] = 0
    if 'shift factor' not in data:
        data['shift factor'] = 1.
    if 'max cyc' not in data:
        data['max cyc'] = 3
    if 'F**2' not in data:
        data['F**2'] = True
        data['minF/sig'] = 0
    if 'Author' not in data:
        data['Author'] = 'no name'
    #end patch

    def SeqSizer():
        
        def OnSelectData(event):
            choices = ['All',]+GetPatternTreeDataNames(G2frame,['PWDR',])
            sel = []
            if 'Seq Data' in data:
                for item in data['Seq Data']:
                    sel.append(choices.index(item))
            names = []
            dlg = wx.MultiChoiceDialog(G2frame,'Select data:','Sequential refinement',choices)
            dlg.SetSelections(sel)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                for i in sel: names.append(choices[i])
                if 'All' in names:
                    names = choices[1:]
                data['Seq Data'] = names                
            dlg.Destroy()
            reverseSel.Enable(True)
            
        def OnReverse(event):
            data['Reverse Seq'] = reverseSel.GetValue()
                    
        seqSizer = wx.BoxSizer(wx.HORIZONTAL)
        seqSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Sequential Refinement Powder Data: '),0,wx.ALIGN_CENTER_VERTICAL)
        selSeqData = wx.Button(G2frame.dataDisplay,-1,label=' Select data')
        selSeqData.Bind(wx.EVT_BUTTON,OnSelectData)
        seqSizer.Add(selSeqData,0,wx.ALIGN_CENTER_VERTICAL)
        seqSizer.Add((5,0),0)
        reverseSel = wx.CheckBox(G2frame.dataDisplay,-1,label=' Reverse order?')
        reverseSel.Bind(wx.EVT_CHECKBOX,OnReverse)
        if 'Seq Data' not in data:
            reverseSel.Enable(False)
        if 'Reverse Seq' in data:
            reverseSel.SetValue(data['Reverse Seq'])
        seqSizer.Add(reverseSel,0,wx.ALIGN_CENTER_VERTICAL)
        return seqSizer
        
    def LSSizer():        
        
        def OnDerivType(event):
            data['deriv type'] = derivSel.GetValue()
            derivSel.SetValue(data['deriv type'])
            wx.CallAfter(UpdateControls,G2frame,data)
            
        def OnConvergence(event):
            try:
                value = max(1.e-9,min(1.0,float(Cnvrg.GetValue())))
            except ValueError:
                value = 0.0001
            data['min dM/M'] = value
            Cnvrg.SetValue('%.2g'%(value))
            
        def OnMaxCycles(event):
            data['max cyc'] = int(maxCyc.GetValue())
            maxCyc.SetValue(str(data['max cyc']))
                        
        def OnFactor(event):
            try:
                value = min(max(float(Factr.GetValue()),0.00001),100.)
            except ValueError:
                value = 1.0
            data['shift factor'] = value
            Factr.SetValue('%.5f'%(value))
            
        def OnFsqRef(event):
            data['F**2'] = fsqRef.GetValue()
        
        def OnMinSig(event):
            try:
                value = min(max(float(minSig.GetValue()),0.),5.)
            except ValueError:
                value = 1.0
            data['minF/sig'] = value
            minSig.SetValue('%.2f'%(value))

        LSSizer = wx.FlexGridSizer(cols=4,vgap=5,hgap=5)
        LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Refinement derivatives: '),0,wx.ALIGN_CENTER_VERTICAL)
        Choice=['analytic Jacobian','numeric','analytic Hessian']
        derivSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['deriv type'],choices=Choice,
            style=wx.CB_READONLY|wx.CB_DROPDOWN)
        derivSel.SetValue(data['deriv type'])
        derivSel.Bind(wx.EVT_COMBOBOX, OnDerivType)
            
        LSSizer.Add(derivSel,0,wx.ALIGN_CENTER_VERTICAL)
        LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Min delta-M/M: '),0,wx.ALIGN_CENTER_VERTICAL)
        Cnvrg = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2g'%(data['min dM/M']),style=wx.TE_PROCESS_ENTER)
        Cnvrg.Bind(wx.EVT_TEXT_ENTER,OnConvergence)
        Cnvrg.Bind(wx.EVT_KILL_FOCUS,OnConvergence)
        LSSizer.Add(Cnvrg,0,wx.ALIGN_CENTER_VERTICAL)
        if 'Hessian' in data['deriv type']:
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Max cycles: '),0,wx.ALIGN_CENTER_VERTICAL)
            Choice = ['0','1','2','3','5','10','15','20']
            maxCyc = wx.ComboBox(parent=G2frame.dataDisplay,value=str(data['max cyc']),choices=Choice,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)
            maxCyc.SetValue(str(data['max cyc']))
            maxCyc.Bind(wx.EVT_COMBOBOX, OnMaxCycles)
            LSSizer.Add(maxCyc,0,wx.ALIGN_CENTER_VERTICAL)
        else:
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Initial shift factor: '),0,wx.ALIGN_CENTER_VERTICAL)
            Factr = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.5f'%(data['shift factor']),style=wx.TE_PROCESS_ENTER)
            Factr.Bind(wx.EVT_TEXT_ENTER,OnFactor)
            Factr.Bind(wx.EVT_KILL_FOCUS,OnFactor)
            LSSizer.Add(Factr,0,wx.ALIGN_CENTER_VERTICAL)
        if G2frame.Sngl:
            LSSizer.Add((1,0),)
            LSSizer.Add((1,0),)
            fsqRef = wx.CheckBox(G2frame.dataDisplay,-1,label='Refine HKLF as F^2? ')
            fsqRef.SetValue(data['F**2'])
            fsqRef.Bind(wx.EVT_CHECKBOX,OnFsqRef)
            LSSizer.Add(fsqRef,0,wx.ALIGN_CENTER_VERTICAL)
            LSSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,label='Min obs/sig (0-5): '),0,wx.ALIGN_CENTER_VERTICAL)
            minSig = wx.TextCtrl(G2frame.dataDisplay,-1,value='%.2f'%(data['minF/sig']),style=wx.TE_PROCESS_ENTER)
            minSig.Bind(wx.EVT_TEXT_ENTER,OnMinSig)
            minSig.Bind(wx.EVT_KILL_FOCUS,OnMinSig)
            LSSizer.Add(minSig,0,wx.ALIGN_CENTER_VERTICAL)
        return LSSizer
        
    def AuthSizer():

        def OnAuthor(event):
            data['Author'] = auth.GetValue()

        Author = data['Author']
        authSizer = wx.BoxSizer(wx.HORIZONTAL)
        authSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' CIF Author (last, first):'),0,wx.ALIGN_CENTER_VERTICAL)
        auth = wx.TextCtrl(G2frame.dataDisplay,-1,value=Author,style=wx.TE_PROCESS_ENTER)
        auth.Bind(wx.EVT_TEXT_ENTER,OnAuthor)
        auth.Bind(wx.EVT_KILL_FOCUS,OnAuthor)
        authSizer.Add(auth,0,wx.ALIGN_CENTER_VERTICAL)
        return authSizer
        
        
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
        Status.SetStatusText('')
    G2frame.dataFrame.SetLabel('Controls')
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    SetDataMenuBar(G2frame,G2frame.dataFrame.ControlsMenu)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0)
    mainSizer.Add(wx.StaticText(G2frame.dataDisplay,label=' Refinement Controls:'),0,wx.ALIGN_CENTER_VERTICAL)    
    mainSizer.Add(LSSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(SeqSizer())
    mainSizer.Add((5,5),0)
    mainSizer.Add(AuthSizer())
    mainSizer.Add((5,5),0)
        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))
     
################################################################################
#####  Comments
################################################################################           
       
def UpdateComments(G2frame,data):                   

    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.dataFrame.SetLabel('Comments')
    G2frame.dataDisplay = wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
        style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_DONTWRAP)
    for line in data:
        G2frame.dataDisplay.AppendText(line+'\n')
    G2frame.dataFrame.setSizePosLeft([400,250])
            
################################################################################
#####  Sequential Results
################################################################################           
       
def UpdateSeqResults(G2frame,data):
    """
    Called when the Sequential Results data tree entry is selected
    to show results from a sequential refinement.
    
    :param wx.Frame G2frame: main GSAS-II data tree windows

    :param dict data: a dictionary containing the following items:  

            * 'histNames' - list of histogram names in order as processed by Sequential Refinement
            * 'varyList' - list of variables - identical over all refinements in sequence
            * 'histName' - dictionaries for all data sets processed, which contains:

              * 'variables'- result[0] from leastsq call
              * 'varyList' - list of variables; same as above
              * 'sig' - esds for variables
              * 'covMatrix' - covariance matrix from individual refinement
              * 'title' - histogram name; same as dict item name
              * 'newAtomDict' - new atom parameters after shifts applied
              * 'newCellDict' - new cell parameters after shifts to A0-A5 applied'
    """
    if not data:
        print 'No sequential refinement results'
        return
    histNames = data['histNames']
       
    def GetSampleParms():
        sampleParmDict = {'Temperature':[],'Pressure':[],'Humidity':[],'Voltage':[],'Force':[],}
        sampleParm = {}
        for name in histNames:
            Id = GetPatternTreeItemId(G2frame,G2frame.root,name)
            sampleData = G2frame.PatternTree.GetItemPyData(GetPatternTreeItemId(G2frame,Id,'Sample Parameters'))
            for item in sampleParmDict:
                sampleParmDict[item].append(sampleData[item])
        for item in sampleParmDict:
            frstValue = sampleParmDict[item][0]
            if np.any(np.array(sampleParmDict[item])-frstValue):
                sampleParm[item] = sampleParmDict[item]            
        return sampleParm
            
    def GetRwps():
        Rwps = []
        for name in histNames:
            Rwps.append(data[name]['Rvals']['Rwp'])
        return Rwps
            
    def GetSigData(parm):
        sigData = []
        for name in histNames:
            sigList = data[name]['sig']
            if colLabels[parm] in atomList:
                sigData.append(sigList[colLabels.index(atomList[colLabels[parm]])-1])
            elif colLabels[parm] in cellList:
                sigData.append(sigList[colLabels.index(cellList[colLabels[parm]])-1])
            else:
                sigData.append(sigList[parm-1])
        return sigData
    
    def Select(event):
        cols = G2frame.dataDisplay.GetSelectedCols()
        rows = G2frame.dataDisplay.GetSelectedRows()
        if cols:
            plotData = []
            plotSig = []
            plotNames = []
            for col in cols:
                plotData.append(G2frame.SeqTable.GetColValues(col))
                if col:     #not Rwp
                    plotSig.append(GetSigData(col))
                else:
                    plotSig.append(0.0)
                plotNames.append(G2frame.SeqTable.GetColLabelValue(col))
            plotData = np.array(plotData)
            G2plt.PlotSeq(G2frame,plotData,plotSig,plotNames,sampleParms)
        elif rows:
            name = histNames[rows[0]]       #only does 1st one selected
            G2plt.PlotCovariance(G2frame,data[name])
            
    def OnSaveSelSeq(event):        
        cols = G2frame.dataDisplay.GetSelectedCols()
        if cols:
            numRows = G2frame.SeqTable.GetNumberRows()
            dataNames = []
            saveNames = [G2frame.SeqTable.GetRowLabelValue(r) for r in range(numRows)]
            saveData = []
            for col in cols:
                dataNames.append(G2frame.SeqTable.GetColLabelValue(col))
                if col:     #not Rwp
                    saveData.append(zip(G2frame.SeqTable.GetColValues(col),GetSigData(col)))
                else:
                    saveData.append(zip(G2frame.SeqTable.GetColValues(col),0.0))
            lenName = len(saveNames[0])
            saveData = np.swapaxes(np.array(saveData),0,1)
            dlg = wx.FileDialog(G2frame, 'Choose text output file for your selection', '.', '', 
                'Text output file (*.txt)|*.txt',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT|wx.CHANGE_DIR)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    SeqTextFile = dlg.GetPath()
                    SeqTextFile = G2IO.FileDlgFixExt(dlg,SeqTextFile)
                    SeqFile = open(SeqTextFile,'w')
                    line = '  %s  '%('name'.center(lenName))
                    for item in dataNames:
                        line += ' %12s %12s '%(item.center(12),'esd'.center(12))
                    line += '\n'
                    SeqFile.write(line)
                    for i,item in enumerate(saveData):
                        line = " '%s' "%(saveNames[i])
                        for val,esd in item:
                            line += ' %12.6f %12.6f '%(val,esd)
                        line += '\n'
                        SeqFile.write(line)
                    SeqFile.close()
            finally:
                dlg.Destroy()
            
               
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    atomList = {}
    newAtomDict = data[histNames[0]]['newAtomDict']
    for item in newAtomDict:
        if item in data['varyList']:
            atomList[newAtomDict[item][0]] = item
    cellList = {}
    newCellDict = data[histNames[0]]['newCellDict']
    for item in newCellDict:
        if item in data['varyList']:
            cellList[newCellDict[item][0]] = item
    sampleParms = GetSampleParms()
    Rwps = GetRwps()
    SetDataMenuBar(G2frame,G2frame.dataFrame.SequentialMenu)
    G2frame.dataFrame.SetLabel('Sequential refinement results')
    G2frame.dataFrame.CreateStatusBar()
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnSaveSelSeq, id=wxID_SAVESEQSEL)
    colLabels = ['Rwp',]+data['varyList']+atomList.keys()+cellList.keys()
    Types = (len(data['varyList']+atomList.keys()+cellList.keys())+1)*[wg.GRID_VALUE_FLOAT,]
    seqList = [[Rwps[i],]+list(data[name]['variables']) for i,name in enumerate(histNames)]
    for i,item in enumerate(seqList):
        newAtomDict = data[histNames[i]]['newAtomDict']
        newCellDict = data[histNames[i]]['newCellDict']
        item += [newAtomDict[atomList[parm]][1] for parm in atomList.keys()]
        item += [newCellDict[cellList[parm]][1] for parm in cellList.keys()]
    G2frame.SeqTable = Table(seqList,colLabels=colLabels,rowLabels=histNames,types=Types)
    G2frame.dataDisplay = GSGrid(parent=G2frame.dataFrame)
    G2frame.dataDisplay.SetTable(G2frame.SeqTable, True)
    G2frame.dataDisplay.EnableEditing(False)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_DCLICK, Select)
    G2frame.dataDisplay.SetRowLabelSize(8*len(histNames[0]))       #pretty arbitrary 8
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(True)
    G2frame.dataFrame.setSizePosLeft([700,350])
        
################################################################################
#####  Main PWDR panel
################################################################################           
       
def UpdatePWHKPlot(G2frame,kind,item):
    '''Called when the histogram main tree entry is called. Displays the
    histogram weight factor, refinement statistics for the histogram
    and the range of data for a simulation.

    Also invokes a plot of the histogram.
    '''
    def onEditSimRange(event):
        'Edit simulation range'
        inp = [
            min(data[1][0]),
            max(data[1][0]),
            None
            ]
        inp[2] = (inp[1] - inp[0])/(len(data[1][0])-1.)
        names = ('start angle', 'end angle', 'step size')
        dictlst = [inp] * len(inp)
        elemlst = range(len(inp))
        dlg = ScrolledMultiEditor(
            G2frame,[inp] * len(inp), range(len(inp)), names,
            header='Edit simulation range',
            minvals=(0.001,0.001,0.0001),
            maxvals=(180.,180.,.1),
            )
        dlg.CenterOnParent()
        val = dlg.ShowModal()
        dlg.Destroy()
        if val != wx.ID_OK: return
        if inp[0] > inp[1]:
            end,start,step = inp
        else:                
            start,end,step = inp
        step = abs(step)
        N = int((end-start)/step)+1
        newdata = np.linspace(start,end,N,True)
        if len(newdata) < 2: return # too small a range - reject
        data[1][0] = newdata
        data[1][1] = np.zeros_like(newdata)
        data[1][2] = np.ones_like(newdata)
        data[1][3] = np.zeros_like(newdata)
        data[1][4] = np.zeros_like(newdata)
        data[1][5] = np.zeros_like(newdata)
        Tmin = newdata[0]
        Tmax = newdata[-1]
        G2frame.PatternTree.SetItemPyData(GetPatternTreeItemId(G2frame,item,'Limits'),
            [(Tmin,Tmax),[Tmin,Tmax]])
        UpdatePWHKPlot(G2frame,kind,item) # redisplay data screen

    def OnErrorAnalysis(event):
        G2plt.PlotDeltSig(G2frame,kind)
        
    def OnWtFactor(event):
        try:
            val = float(wtval.GetValue())
        except ValueError:
            val = data[0]['wtFactor']
        data[0]['wtFactor'] = val
        wtval.SetValue('%.3f'%(val))
           
    data = G2frame.PatternTree.GetItemPyData(item)
    if 'wtFactor' not in data[0]:
        data[0] = {'wtFactor':1.0}
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    SetDataMenuBar(G2frame,G2frame.dataFrame.ErrorMenu)
    G2frame.dataFrame.Bind(wx.EVT_MENU,OnErrorAnalysis, id=wxID_PWDANALYSIS)
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),)
    wtSizer = wx.BoxSizer(wx.HORIZONTAL)
    wtSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' Weight factor: '),0,wx.ALIGN_CENTER_VERTICAL)
    wtval = wx.TextCtrl(G2frame.dataDisplay,-1,'%.3f'%(data[0]['wtFactor']),style=wx.TE_PROCESS_ENTER)
    wtval.Bind(wx.EVT_TEXT_ENTER,OnWtFactor)
    wtval.Bind(wx.EVT_KILL_FOCUS,OnWtFactor)
    wtSizer.Add(wtval,0,wx.ALIGN_CENTER_VERTICAL)
    mainSizer.Add(wtSizer)
    if data[0].get('Dummy'):
        simSizer = wx.BoxSizer(wx.HORIZONTAL)
        Tmin = min(data[1][0])
        Tmax = max(data[1][0])
        num = len(data[1][0])
        step = (Tmax - Tmin)/(num-1)
        t = u'2\u03b8' # 2theta
        lbl =  u'Simulation range: {:.2f} to {:.2f} {:s}\nwith {:.4f} steps ({:d} points)'
        lbl += u'\n(Edit range resets observed intensities).'
        lbl = lbl.format(Tmin,Tmax,t,step,num)
        simSizer.Add(wx.StaticText(G2frame.dataDisplay,wx.ID_ANY,lbl),
                    0,wx.ALIGN_CENTER_VERTICAL)
        but = wx.Button(G2frame.dataDisplay,wx.ID_ANY,"Edit range")
        but.Bind(wx.EVT_BUTTON,onEditSimRange)
        simSizer.Add(but,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(simSizer)
    if 'Nobs' in data[0]:
        mainSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,
            ' Data residual wR: %.3f%% on %d observations'%(data[0]['wR'],data[0]['Nobs'])))
        for value in data[0]:
            if 'Nref' in value:
                mainSizer.Add((5,5),)
                pfx = value.split('Nref')[0]
                name = data[0][pfx.split(':')[0]+'::Name']
                mainSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,' For phase '+name+':'))
                mainSizer.Add(wx.StaticText(G2frame.dataDisplay,-1,
                    ' Unweighted phase residuals RF^2: %.3f%%, RF: %.3f%% on %d reflections  '% \
                    (data[0][pfx+'Rf^2'],data[0][pfx+'Rf'],data[0][value])))
    mainSizer.Add((5,5),)
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    Size = mainSizer.Fit(G2frame.dataFrame)
    Size[1] += 10
    G2frame.dataFrame.setSizePosLeft(Size)
    G2frame.PatternTree.SetItemPyData(item,data)
    if kind == 'PWDR':
        G2plt.PlotPatterns(G2frame,newPlot=True)
    elif kind == 'HKLF':
        G2plt.PlotSngl(G2frame,newPlot=True)
                 
################################################################################
#####  HKLF controls
################################################################################           
       
def UpdateHKLControls(G2frame,data):
    '''Needs a doc string
    '''
    
    def OnScaleSlider(event):
        scale = int(scaleSel.GetValue())/1000.
        scaleSel.SetValue(int(scale*1000.))
        data['Scale'] = scale*1.
        G2plt.PlotSngl(G2frame)
        
    def OnLayerSlider(event):
        layer = layerSel.GetValue()
        data['Layer'] = layer
        G2plt.PlotSngl(G2frame)
        
    def OnSelZone(event):
        data['Zone'] = zoneSel.GetValue()
        izone = zones.index(data['Zone'])
        layerSel.SetRange(maxValue=HKLmax[izone],minValue=HKLmin[izone])
        G2plt.PlotSngl(G2frame,newPlot=True)
        
    def OnSelType(event):
        data['Type'] = typeSel.GetValue()
        G2plt.PlotSngl(G2frame)
        
    def SetStatusLine():
        Status.SetStatusText("")
                                      
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    SetStatusLine()
    zones = ['100','010','001']
    HKLmax = data['HKLmax']
    HKLmin = data['HKLmin']
    typeChoices = ['Fosq','Fo','|DFsq|/sig','|DFsq|>sig','|DFsq|>3sig']
    G2frame.dataDisplay = wx.Panel(G2frame.dataFrame)
    SetDataMenuBar(G2frame)
    G2frame.dataFrame.SetTitle('HKL Plot Controls')
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,10),0)
    
    scaleSizer = wx.BoxSizer(wx.HORIZONTAL)
    scaleSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Scale'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    scaleSel = wx.Slider(parent=G2frame.dataDisplay,maxValue=1000,minValue=1,
        style=wx.SL_HORIZONTAL,value=int(data['Scale']*10))
    scaleSizer.Add(scaleSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    scaleSel.SetLineSize(10)
    scaleSel.SetPageSize(10)
    scaleSel.Bind(wx.EVT_SLIDER, OnScaleSlider)
    mainSizer.Add(scaleSizer,0,wx.EXPAND|wx.RIGHT)
    mainSizer.Add((0,10),0)    
    
    zoneSizer = wx.BoxSizer(wx.HORIZONTAL)
    zoneSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Zone  '),0,
        wx.ALIGN_CENTER_VERTICAL)
    zoneSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['Zone'],choices=['100','010','001'],
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    zoneSel.Bind(wx.EVT_COMBOBOX, OnSelZone)
    zoneSizer.Add(zoneSel,0,wx.ALIGN_CENTER_VERTICAL)
    zoneSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Plot type  '),0,
        wx.ALIGN_CENTER_VERTICAL)        
    typeSel = wx.ComboBox(parent=G2frame.dataDisplay,value=data['Type'],choices=typeChoices,
        style=wx.CB_READONLY|wx.CB_DROPDOWN)
    typeSel.Bind(wx.EVT_COMBOBOX, OnSelType)
    zoneSizer.Add(typeSel,0,wx.ALIGN_CENTER_VERTICAL)
    zoneSizer.Add((10,0),0)    
    mainSizer.Add(zoneSizer,0,wx.EXPAND|wx.RIGHT)
    mainSizer.Add((0,10),0)    
        
    izone = zones.index(data['Zone'])
    layerSizer = wx.BoxSizer(wx.HORIZONTAL)
    layerSizer.Add(wx.StaticText(parent=G2frame.dataDisplay,label=' Layer'),0,
        wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
    layerSel = wx.Slider(parent=G2frame.dataDisplay,maxValue=HKLmax[izone],minValue=HKLmin[izone],
        style=wx.SL_HORIZONTAL|wx.SL_AUTOTICKS|wx.SL_LABELS,value=0)
    layerSel.SetLineSize(1)
    layerSel.SetPageSize(1)
    layerSel.Bind(wx.EVT_SLIDER, OnLayerSlider)    
    layerSizer.Add(layerSel,1,wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
    layerSizer.Add((10,0),0)    
    mainSizer.Add(layerSizer,1,wx.EXPAND|wx.RIGHT)

        
    mainSizer.Layout()    
    G2frame.dataDisplay.SetSizer(mainSizer)
    G2frame.dataDisplay.SetSize(mainSizer.Fit(G2frame.dataFrame))
    G2frame.dataFrame.setSizePosLeft(mainSizer.Fit(G2frame.dataFrame))

################################################################################
#####  Pattern tree routines
################################################################################           
       
def GetPatternTreeDataNames(G2frame,dataTypes):
    '''Needs a doc string
    '''
    names = []
    item, cookie = G2frame.PatternTree.GetFirstChild(G2frame.root)        
    while item:
        name = G2frame.PatternTree.GetItemText(item)
        if name[:4] in dataTypes:
            names.append(name)
        item, cookie = G2frame.PatternTree.GetNextChild(G2frame.root, cookie)
    return names
                          
def GetPatternTreeItemId(G2frame, parentId, itemText):
    '''Needs a doc string
    '''
    item, cookie = G2frame.PatternTree.GetFirstChild(parentId)
    while item:
        if G2frame.PatternTree.GetItemText(item) == itemText:
            return item
        item, cookie = G2frame.PatternTree.GetNextChild(parentId, cookie)
    return 0                

def MovePatternTreeToGrid(G2frame,item):
    '''Needs a doc string
    '''
    
#    print G2frame.PatternTree.GetItemText(item)
    
    oldPage = None # will be set later if already on a Phase item
    if G2frame.dataFrame:
        SetDataMenuBar(G2frame)
        if G2frame.dataFrame.GetLabel() == 'Comments':
            try:
                data = [G2frame.dataDisplay.GetValue()]
                G2frame.dataDisplay.Clear() 
                Id = GetPatternTreeItemId(G2frame,G2frame.root, 'Comments')
                if Id: G2frame.PatternTree.SetItemPyData(Id,data)
            except:     #clumsy but avoids dead window problem when opening another project
                pass
        elif G2frame.dataFrame.GetLabel() == 'Notebook':
            try:
                data = [G2frame.dataDisplay.GetValue()]
                G2frame.dataDisplay.Clear() 
                Id = GetPatternTreeItemId(G2frame,G2frame.root, 'Notebook')
                if Id: G2frame.PatternTree.SetItemPyData(Id,data)
            except:     #clumsy but avoids dead window problem when opening another project
                pass
        elif 'Phase Data for' in G2frame.dataFrame.GetLabel():
            if G2frame.dataDisplay: 
                oldPage = G2frame.dataDisplay.GetSelection()
        G2frame.dataFrame.Clear()
        G2frame.dataFrame.SetLabel('')
    else:
        #create the frame for the data item window
        G2frame.dataFrame = DataFrame(parent=G2frame.mainPanel,frame=G2frame)
        G2frame.dataFrame.PhaseUserSize = None
        
    G2frame.dataFrame.Raise()            
    G2frame.PickId = 0
    parentID = G2frame.root
    for i in G2frame.ExportPattern: i.Enable(False)
    defWid = [250,150]
    if item != G2frame.root:
        parentID = G2frame.PatternTree.GetItemParent(item)
    if G2frame.PatternTree.GetItemParent(item) == G2frame.root:
        G2frame.PatternId = item
        G2frame.PickId = item
        if G2frame.PatternTree.GetItemText(item) == 'Notebook':
            SetDataMenuBar(G2frame,G2frame.dataFrame.DataNotebookMenu)
            G2frame.PatternId = 0
            for i in G2frame.ExportPattern: i.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateNotebook(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Controls':
            G2frame.PatternId = 0
            for i in G2frame.ExportPattern: i.Enable(False)
            data = G2frame.PatternTree.GetItemPyData(item)
            if not data:           #fill in defaults
                data = {
                    #least squares controls
                    'deriv type':'analytic Hessian','min dM/M':0.0001,'shift factor':1.0,'max cyc':3}
                G2frame.PatternTree.SetItemPyData(item,data)                             
            for i in G2frame.Refine: i.Enable(True)
            for i in G2frame.SeqRefine: i.Enable(True)
            UpdateControls(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Sequential results':
            data = G2frame.PatternTree.GetItemPyData(item)
            UpdateSeqResults(G2frame,data)            
        elif G2frame.PatternTree.GetItemText(item) == 'Covariance':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2frame.dataFrame.setSizePosLeft(defWid)
            text = ''
            if 'Rvals' in data:
                Nvars = len(data['varyList'])
                Rvals = data['Rvals']
                text = '\nFinal residuals: \nwR = %.3f%% \nchi**2 = %.1f \nGOF = %.2f'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF'])
                text += '\nNobs = %d \nNvals = %d'%(Rvals['Nobs'],Nvars)
                if 'lamMax' in Rvals:
                    text += '\nlog10 MaxLambda = %.1f'%(np.log10(Rvals['lamMax']))
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                value='See plot window for covariance display'+text,style=wx.TE_MULTILINE)
            G2plt.PlotCovariance(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Constraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2cnstG.UpdateConstraints(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Rigid bodies':
            data = G2frame.PatternTree.GetItemPyData(item)
            G2cnstG.UpdateRigidBodies(G2frame,data)
        elif G2frame.PatternTree.GetItemText(item) == 'Restraints':
            data = G2frame.PatternTree.GetItemPyData(item)
            Phases = G2frame.GetPhaseData()
            phase = ''
            phaseName = ''
            if Phases:
                phaseName = Phases.keys()[0]
            G2frame.dataFrame.setSizePosLeft(defWid)
            G2restG.UpdateRestraints(G2frame,data,Phases,phaseName)
        elif 'IMG' in G2frame.PatternTree.GetItemText(item):
            G2frame.Image = item
            G2plt.PlotImage(G2frame,newPlot=True)
        elif 'PKS' in G2frame.PatternTree.GetItemText(item):
            G2plt.PlotPowderLines(G2frame)
        elif 'PWDR' in G2frame.PatternTree.GetItemText(item):
            for i in G2frame.ExportPattern: i.Enable(True)
            UpdatePWHKPlot(G2frame,'PWDR',item)
        elif 'HKLF' in G2frame.PatternTree.GetItemText(item):
            G2frame.Sngl = item
            UpdatePWHKPlot(G2frame,'HKLF',item)
        elif 'PDF' in G2frame.PatternTree.GetItemText(item):
            G2frame.PatternId = item
            for i in G2frame.ExportPDF: i.Enable(True)
            G2plt.PlotISFG(G2frame,type='S(Q)')
        elif G2frame.PatternTree.GetItemText(item) == 'Phases':
            G2frame.dataFrame.setSizePosLeft(defWid)
            wx.TextCtrl(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize(),
                value='Select one phase to see its parameters')            
    elif 'I(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='I(Q)',newPlot=True)
    elif 'S(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='S(Q)',newPlot=True)
    elif 'F(Q)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='F(Q)',newPlot=True)
    elif 'G(R)' in G2frame.PatternTree.GetItemText(item):
        G2frame.PickId = item
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2plt.PlotISFG(G2frame,type='G(R)',newPlot=True)            
    elif G2frame.PatternTree.GetItemText(parentID) == 'Phases':
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2phG.UpdatePhaseData(G2frame,item,data,oldPage)
    elif G2frame.PatternTree.GetItemText(item) == 'Comments':
        SetDataMenuBar(G2frame,G2frame.dataFrame.DataCommentsMenu)
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        UpdateComments(G2frame,data)
    elif G2frame.PatternTree.GetItemText(item) == 'Image Controls':
        G2frame.dataFrame.SetTitle('Image Controls')
        G2frame.PickId = item
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        masks = G2frame.PatternTree.GetItemPyData(
            GetPatternTreeItemId(G2frame,G2frame.Image, 'Masks'))
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateImageControls(G2frame,data,masks)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Masks':
        G2frame.dataFrame.SetTitle('Masks')
        G2frame.PickId = item
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateMasks(G2frame,data)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Stress/Strain':
        G2frame.dataFrame.SetTitle('Stress/Strain')
        G2frame.PickId = item
        G2frame.Image = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        G2imG.UpdateStressStrain(G2frame,data)
        G2plt.PlotImage(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'HKL Plot Controls':
        G2frame.PickId = item
        G2frame.Sngl = G2frame.PatternTree.GetItemParent(item)
        data = G2frame.PatternTree.GetItemPyData(item)
        UpdateHKLControls(G2frame,data)
        G2plt.PlotSngl(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'PDF Controls':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPDF: i.Enable(True)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePDFGrid(G2frame,data)
        G2plt.PlotISFG(G2frame,type='I(Q)')
        G2plt.PlotISFG(G2frame,type='S(Q)')
        G2plt.PlotISFG(G2frame,type='F(Q)')
        G2plt.PlotISFG(G2frame,type='G(R)')
    elif G2frame.PatternTree.GetItemText(item) == 'Peak List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdatePeakGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Background':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateBackground(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Limits':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateLimitsGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Instrument Parameters':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)[0]
        G2pdG.UpdateInstrumentGrid(G2frame,data)
        G2plt.PlotPeakWidths(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Sample Parameters':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)

        if 'Temperature' not in data:           #temp fix for old gpx files
            data = {'Scale':[1.0,True],'Type':'Debye-Scherrer','Absorption':[0.0,False],'DisplaceX':[0.0,False],
                'DisplaceY':[0.0,False],'Diffuse':[],'Temperature':300.,'Pressure':1.0,'Humidity':0.0,'Voltage':0.0,
                'Force':0.0,'Gonio. radius':200.0}
            G2frame.PatternTree.SetItemPyData(item,data)
    
        G2pdG.UpdateSampleGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Index Peak List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        for i in G2frame.ExportPeakList: i.Enable(True)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2pdG.UpdateIndexPeaksGrid(G2frame,data)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Unit Cells List':
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        if not data:
            data.append([0,0.0,4,25.0,0,'P1',1,1,1,90,90,90]) #zero error flag, zero value, max Nc/No, start volume
            data.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0])      #Bravais lattice flags
            data.append([])                                 #empty cell list
            data.append([])                                 #empty dmin
            G2frame.PatternTree.SetItemPyData(item,data)                             
        G2pdG.UpdateUnitCellsGrid(G2frame,data)
        if 'PKS' in G2frame.PatternTree.GetItemText(G2frame.PatternId):
            G2plt.PlotPowderLines(G2frame)
        else:
            G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Reflection Lists':   #powder reflections
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        G2frame.PickId = item
        data = G2frame.PatternTree.GetItemPyData(item)
        G2frame.RefList = ''
        if len(data):
            G2frame.RefList = data.keys()[0]
        G2pdG.UpdateReflectionGrid(G2frame,data)
        G2plt.PlotPatterns(G2frame)
    elif G2frame.PatternTree.GetItemText(item) == 'Reflection List':    #HKLF reflections
        G2frame.PatternId = G2frame.PatternTree.GetItemParent(item)
        name = G2frame.PatternTree.GetItemText(G2frame.PatternId)
        data = G2frame.PatternTree.GetItemPyData(G2frame.PatternId)
        G2pdG.UpdateReflectionGrid(G2frame,data,HKLF=True,Name=name)

def SetDataMenuBar(G2frame,menu=None):
    '''Set the menu for the data frame. On the Mac put this
    menu for the data tree window instead.

    Note that data frame items do not have menus, for these (menu=None)
    display a blank menu or on the Mac display the standard menu for
    the data tree window.
    '''
    if sys.platform == "darwin":
        if menu is None:
            G2frame.SetMenuBar(G2frame.GSASIIMenu)
        else:
            G2frame.SetMenuBar(menu)
    else:
        if menu is None:
            G2frame.dataFrame.SetMenuBar(G2frame.dataFrame.BlankMenu)
        else:
            G2frame.dataFrame.SetMenuBar(menu)

def HorizontalLine(sizer,parent):
    '''Draws a horizontal line as wide as the window.
    This shows up on the Mac as a very thin line, no matter what I do
    '''
    line = wx.StaticLine(parent,-1, size=(-1,3), style=wx.LI_HORIZONTAL)
    sizer.Add(line, 0, wx.EXPAND|wx.ALIGN_CENTER|wx.ALL, 10)

if __name__ == '__main__':
    # test ScrolledMultiEditor
    app = wx.PySimpleApp()
    frm = wx.Frame(None) # create a frame
    frm.Show(True)
    Data1 = {
        'Order':1,
        'omega':'string',
        'chi':2.0,
        'phi':'',
        }
    elemlst = sorted(Data1.keys())
    postlbl = sorted(Data1.keys())
    dictlst = len(elemlst)*[Data1,]

    Data2 = list(range(100))
    elemlst += range(2,6)
    postlbl += range(2,6)
    dictlst += len(range(2,6))*[Data2,]

    prelbl = range(len(elemlst))
    postlbl[1] = "a very long label for the 2nd item to force a horiz. scrollbar"
    header="""This is a longer\nmultiline and perhaps silly header"""
    dlg = ScrolledMultiEditor(frm,dictlst,elemlst,prelbl,postlbl,
                              header=header,CopyButton=True)
    print Data1
    if dlg.ShowModal() == wx.ID_OK:
        for d,k in zip(dictlst,elemlst):
            print k,d[k]
    dlg.Destroy()
    if CallScrolledMultiEditor(frm,dictlst,elemlst,prelbl,postlbl,
                               header=header):
        for d,k in zip(dictlst,elemlst):
            print k,d[k]

#app.MainLoop()
