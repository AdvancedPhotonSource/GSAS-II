# -*- coding: utf-8 -*-
#GSASIIexprGUI - Expression Definition and Evaluation
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIexprGUI: Expression Handling*
-------------------------------------

This module defines a class for defining an expression in terms of values
in a parameter dictionary via a wx.Dialog. The dialog creates a :class:`GSASII.ExpressionObj`
which is used to evaluate the expression against a supplied parameter dictionary.

The expression is parsed to find variables used in the expression and then
the user is asked to assign parameters from the dictionary to each variable.

Default expressions are read from file DefaultExpressions.txt using
:func:`GSASIIpath.LoadConfigFile`.

'''
import re
import sys
import wx
import os.path
import wx.lib.scrolledpanel as wxscroll
import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIgrid as G2gd
import GSASIIctrls as G2G
import GSASIIpy3 as G2py3
import GSASIIobj as G2obj

def IndexParmDict(parmDict,wildcard):
    '''Separate the parameters in parmDict into list of keys by parameter
    type.
    
    :param dict parmDict: a dict with GSAS-II parameters
    :param bool wildcard: True if wildcard versions of parameters should
      be generated and added to the lists
    :returns: a dict of lists where key 1 is a list of phase parameters,
      2 is histogram/phase parms, 3 is histogram parms and 4 are global parameters
    '''
    prex = re.compile('[0-9]+::.*')
    hrex = re.compile(':[0-9]+:.*')
    parmLists = {}
    for i in (1,2,3,4):
        parmLists[i] = []
    for i in sorted(parmDict.keys()):
        if i.startswith("::") or i.find(':') == -1: # globals
            parmLists[4].append(i)
        elif prex.match(i):
            parmLists[1].append(i)
        elif hrex.match(i):
            parmLists[3].append(i)
        else:
            parmLists[2].append(i)
    if wildcard:
        for i in (1,2,3,4):
            parmLists[i] += G2obj.GenWildCard(parmLists[i]) # generate wildcard versions
    for i in (1,2,3,4):
        parmLists[i].sort()
    return parmLists

#==========================================================================
defaultExpressions = None
def LoadDefaultExpressions():
    '''Read a configuration file with default expressions from all files named
    DefaultExpressions.txt found in the path. Duplicates are removed and
    expressions are sorted alphabetically
    '''
    global defaultExpressions
    if defaultExpressions is not None: return # run this routine only once
    defaultExpressions = sorted(list(set(GSASIIpath.LoadConfigFile('DefaultExpressions.txt'))))
    
#==========================================================================
class ExpressionDialog(wx.Dialog):
    '''A wx.Dialog that allows a user to input an arbitrary expression
    to be evaluated and possibly minimized.

    To do this, the user assigns a new (free) or existing
    GSAS-II parameter to each parameter label used in the expression.
    The free parameters can optionally be designated to be refined. 
    For example, is an expression is used such as::

    'A*np.exp(-B/C)'

    then A, B and C can each be assigned as Free parameter with a user-selected
    value or to any existing GSAS-II variable in the parameter dictionary.
    As the expression is entered it is checked for validity.

    After the :class:`ExpressionDialog` object is created, use :meth:`Show` to
    run it and obtain a :class:`GSASIIobj.ExpressionObj` object with the user
    input.

    :param wx.Frame parent: The parent of the Dialog. Can be None,
      but better is to provide the name of the Frame where the dialog
      is called.
    :param dict parmDict: a dict with defined parameters and their values. Each value
      may be a list with parameter values and a refine flag or may just contain
      the parameter value (non-float/int values in dict are ignored)
    :param str exprObj: a :class:`GSASIIobj.ExpressionObj` object with an expression and
      label assignments or None (default)
    :param str wintitle: String placed on title bar of dialog; 
      defaults to "Expression Editor"
    :param str header: String placed at top of dialog to tell the user
      what they will do here; default is "Enter restraint expression here"
    :param bool fit: determines if the expression will be used in fitting (default=True).
      If set to False, and refinement flags are not shown
      and Free parameters are not offered as an assignment option.
    :param str VarLabel: an optional variable label to include before the expression
      input. Ignored if None (default)
    :param list depVarDict: a dict of choices for the dependent variable to be
      fitted to the expression and their values. Ignored if None (default).
    :param list ExtraButton: a list with two terms that define [0]: the label
      for an extra button and [1] the callback routine to be used when the
      button is pressed. The button will only be enabled when the OK button can be
      used (meaning the equation/expression is valid). The default is None, meaning this
      will not be used.
    :param list usedVars: defines a list of previously used variable names. These names
      will not be reused as defaults for new free variables.
      (The default is an empty list). 
    '''
    def __init__(self, parent, parmDict, exprObj=None,
                 header='Enter restraint expression here',
                 wintitle='Expression Editor',
                 fit=True,VarLabel=None,depVarDict=None,
                 ExtraButton=None,usedVars=[]):
        self.fit = fit
        self.depVarDict = depVarDict
        'dict for dependent variables'
        self.parmDict = {}
        '''A copy of the G2 parameter dict (parmDict) except only numerical
        values are included and only the value (not the vary flag, if present)
        is included.
        '''
        self.exprVarLst = []
        '''A list containing the variables utilized in the current expression.
        Placed into a :class:`GSASIIobj.ExpressionObj` object when the dialog is closed
        with "OK", saving any changes.
        '''
        self.varSelect = {}
        '''A dict that shows the variable type for each label
        found in the expression.

        * If the value is None or is not defined, the value is not assigned.
        * If the value is 0, then the varible is "free" -- a new refineable
          parameter.
        * Values above 1 determine what variables will be shown
          when the option is selected. 
        '''
        self.varName = {}
        'Name assigned to each variable'
        self.varValue = {}
        'Value for a variable (Free parameters only)'
        self.varRefflag = {}
        'Refinement flag for a variable (Free parameters only)'
        self.expr = ''
        'Expression as a text string'
        self.dependentVar = None
        'name for dependent variable selection, when depVarDict is specified'
        self.usedVars = usedVars
        'variable names that have been used and should not be reused by default'
        defSize = (620,340) # seems like a good size
        'Starting size for dialog'

        # process dictionary of values and create an index
        for key in parmDict:
            try: # deal with values that are in lists
                val = parmDict[key][0]
            except KeyError:
                continue # there were dicts in parmDict (should be gone now)
            except (TypeError,IndexError):
                val = parmDict[key]
            if isinstance(val, basestring): continue
            try:
                self.parmDict[key] = float(val)
            except:
                pass
        # separate the variables by type
        self.parmLists = IndexParmDict(self.parmDict,self.fit)
        self.timer = wx.Timer()
        self.timer.Bind(wx.EVT_TIMER,self.OnValidate)
        LoadDefaultExpressions()
        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, parent, wx.ID_ANY, wintitle, style=style, size=defSize)
        self.mainsizer = wx.BoxSizer(wx.VERTICAL)
        label = wx.StaticText(self,  wx.ID_ANY, header)
        self.mainsizer.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)

        self.exsizer = wx.BoxSizer(wx.HORIZONTAL)
        if VarLabel:
            label = wx.StaticText(self,  wx.ID_ANY, VarLabel + ' = ')
            self.exsizer.Add(label, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 0)
        elif depVarDict:
            self.depParmLists = IndexParmDict(self.depVarDict,False)
            choices = ['','Phase','Hist./Phase','Hist.','Global']
            for i in range(1,len(choices)): # remove empty menus from choice list
                if not len(self.depParmLists[i]): choices[i] = ''
            self.depChoices = [i for i in range(len(choices)) if choices[i]]
            choice = wx.Choice(
                self, wx.ID_ANY,
                choices = [choices[i] for i in self.depChoices]
                )
            choice.SetSelection(wx.NOT_FOUND)
            choice.Bind(wx.EVT_CHOICE,self.OnDepChoice)
            self.exsizer.Add(choice, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 0)
            self.exsizer.Add((5,5))
            self.depLabel = wx.StaticText(self,  wx.ID_ANY, ' ')
            self.exsizer.Add(self.depLabel, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 0)
            label = wx.StaticText(self,  wx.ID_ANY, ' = ')
            self.exsizer.Add(label, 0, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 0)

        #self.exCtrl = wx.TextCtrl(self,  wx.ID_ANY, size=(150,-1),style=wx.TE_PROCESS_ENTER)
        self.exCtrl = wx.ComboBox(self, wx.ID_ANY, "", (90, 50), (160, -1),
                                  defaultExpressions,
                                  wx.CB_DROPDOWN| wx.TE_PROCESS_ENTER
                         )
        self.exCtrl.Bind(wx.EVT_CHAR, self.OnChar)
        self.exCtrl.Bind(wx.EVT_COMBOBOX, self.OnValidate)
        self.exCtrl.Bind(wx.EVT_TEXT_ENTER, self.OnValidate)
        self.exsizer.Add(self.exCtrl, 1, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 0)
        #self.mainsizer.Add(self.exCtrl, 0, wx.ALL|wx.EXPAND, 5)
        self.mainsizer.Add(self.exsizer, 0, wx.ALL|wx.EXPAND, 5)
        self.mainsizer.Add((-1,5),0,wx.EXPAND,1)

        evalSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.mainsizer.Add(evalSizer,0,wx.ALL|wx.EXPAND,0)
        btn = wx.Button(self, wx.ID_ANY,"Validate")
        btn.Bind(wx.EVT_BUTTON,self.OnValidate)
        evalSizer.Add(btn,0,wx.LEFT|wx.RIGHT,5)
        self.result = wx.StaticText(self,  wx.ID_ANY, '')
        evalSizer.Add(self.result, 0, wx.ALIGN_CENTRE|wx.ALL, 5)

        self.varSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.mainsizer.Add(self.varSizer,1,wx.ALL|wx.EXPAND,1)
        self.mainsizer.Add((-1,5),0,wx.EXPAND,1)

        bSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnsizer = wx.StdDialogButtonSizer()
        if ExtraButton:
            self.ExtraBtn = wx.Button(self, wx.ID_ANY, ExtraButton[0])
            self.ExtraBtn.Bind(wx.EVT_BUTTON,self.OnExtra)
            self.ExtraCallBack = ExtraButton[1]
            self.ExtraBtn.Disable()
            bSizer.Add(self.ExtraBtn, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 2)
        else:
            self.ExtraBtn = None
        bSizer.Add((1,1), 1, wx.ALIGN_CENTER|wx.ALL|wx.EXPAND, 0)
        self.OKbtn = wx.Button(self, wx.ID_OK)
        self.OKbtn.SetDefault()
        self.OKbtn.Disable()
        btnsizer.AddButton(self.OKbtn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        bSizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        self.mainsizer.Add(bSizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL|wx.EXPAND, 5)
        self.SetSizer(self.mainsizer)
        self.CenterOnParent()
        if exprObj:
            self.expr = exprObj.EditExpression(
                self.exprVarLst,
                self.varSelect,
                self.varName,
                self.varValue,
                self.varRefflag,
                )
            # set the initial value for the dependent value
            if self.depVarDict:
                var = exprObj.GetDepVar()
                if var in self.depVarDict:
                    self.depLabel.SetLabel(var)
                    self.dependentVar = var
                    
        self.exCtrl.SetValue(self.expr)
        self.OnValidate(None)
        self.SetMinSize(defSize) 
        #self.errbox.SetAutoLayout(1)
        #self.errbox.SetupScrolling()
        #self.varbox.SetAutoLayout(1)
        #self.varbox.SetupScrolling()
        #self.mainsizer.Fit(self)

    def OnExtra(self,event):
        exprObj = G2obj.ExpressionObj()
        exprObj.LoadExpression(
            self.expr,
            self.exprVarLst,
            self.varSelect,
            self.varName,
            self.varValue,
            self.varRefflag,
            )
        if self.depVarDict:
            exprObj.SetDepVar(self.dependentVar)
        self.ExtraCallBack(exprObj)
        # put results back into displayed dialog
        resDict = dict(exprObj.GetVariedVarVal())
        for v,var in self.varName.items():
            varname = "::" + var.lstrip(':').replace(' ','_').replace(':',';')
            val =  resDict.get(varname)
            if val:
                self.varValue[v] = val
        wx.CallAfter(self.ShowVars)

    def Show(self,mode=True):
        '''Call to use the dialog after it is created.

        :returns: None (On Cancel) or a new :class:`~GSASIIobj.ExpressionObj`
        '''
        self.Layout()
        #self.mainsizer.Fit(self)
        self.SendSizeEvent() # force repaint
        if self.ShowModal() == wx.ID_OK:
            # store the edit results in the object and return it
            exprObj = G2obj.ExpressionObj()
            exprObj.LoadExpression(
                self.expr,
                self.exprVarLst,
                self.varSelect,
                self.varName,
                self.varValue,
                self.varRefflag,
                )
            if self.depVarDict:
                exprObj.SetDepVar(self.dependentVar)
            return exprObj
        else:
            return None
        
    def setEvalResult(self,msg):
        'Show a string in the expression result area'
        self.result.SetLabel(msg)

    def RestartTimer(self):
        '''Cancels any running timer and starts a new one.
        The timer causes a check of syntax after 2 seconds unless there is further input.
        Disables the OK button until a validity check is complete.
        '''
        if self.timer.IsRunning():
            self.timer.Stop()
        self.timer.Start(2000,oneShot=True)
        
    def OnChar(self,event):
        '''Called as each character is entered. Cancels any running timer
        and starts a new one. The timer causes a check of syntax after 2 seconds
        without input.
        Disables the OK button until a validity check is complete.
        '''
        self.RestartTimer()
        self.OKbtn.Disable()
        if self.ExtraBtn: self.ExtraBtn.Disable()
        event.Skip()
        return
    
    def CheckVars(self):
        '''Check that appropriate variables are defined for each
        symbol used in :data:`self.expr`

        :returns: a text error message or None if all needed input is present       
        '''
        invalid = 0
        for v in self.exprVarLst:
            if self.varSelect.get(v) is None:
                invalid += 1
        if invalid==1:
            return '(a variable is not assigned)'
        elif invalid:
            return '('+str(invalid)+' variables are not assigned)'
        msg = ''
        for v in self.exprVarLst:
            varname = self.varName.get(v)
            if not varname:
                invalid += 1
                if msg: msg += "; "
                msg += 'No variable for '+str(v)
            elif self.varSelect.get(v) > 0:
               if '*' in varname:
                   l = G2obj.LookupWildCard(varname,self.parmDict.keys())
                   if len(l) == 0:
                       invalid += 1
                       if msg: msg += "; "
                       msg += 'No variables match '+str(varname)
               elif varname not in self.parmDict.keys():
                   invalid += 1
                   if msg: msg += "; "
                   msg += 'No variables match '+str(varname)
            else:
                # value assignment: this check is likely unneeded
                val = self.varValue.get(v)
                try:
                    float(val)
                except ValueError,TypeError:
                    invalid += 1
                    if msg: msg += "; "
                    if val is None:
                        msg += 'No value for '+str(v)
                    else:
                        msg += 'Value '+str(val)+' invalid for '+str(v)
        if invalid:
            return '('+msg+')'        
        return

    def ShowVars(self):
        # create widgets to associate vars with labels and/or show messages
        self.varSizer.Clear(True)
        self.errbox = wxscroll.ScrolledPanel(self,style=wx.HSCROLL)
        self.errbox.SetMinSize((100,130))
        self.varSizer.Add(self.errbox,0,wx.ALL|wx.EXPAND,1)
        self.varbox = wxscroll.ScrolledPanel(self,style=wx.HSCROLL)
        self.varSizer.Add(self.varbox,1,wx.ALL|wx.EXPAND,1)
        Siz = wx.BoxSizer(wx.VERTICAL)
        Siz.Add(
            wx.StaticText(self.varbox,wx.ID_ANY,
                          'Assignment of variables to labels:'),
            0,wx.EXPAND|wx.ALIGN_CENTER,0)
        GridSiz = wx.FlexGridSizer(0,5,2,2)
        GridSiz.Add(
            wx.StaticText(self.varbox,wx.ID_ANY,'label',style=wx.CENTER),
            0,wx.ALIGN_CENTER)
        lbls = ('varib. type\nselection','variable\nname','value')
        choices = ['Free','Phase','Hist./Phase','Hist.','Global']
        if self.fit:
            lbls += ('refine\nflag',)
        else:
            lbls += ('',)
            choices[0] = ''
        for i in range(1,len(choices)): # remove empty menus from choice list
            if not len(self.parmLists[i]): choices[i] = ''
        self.AllowedChoices = [i for i in range(len(choices)) if choices[i]]
        for lbl in lbls:
            w = wx.StaticText(self.varbox,wx.ID_ANY,lbl,style=wx.CENTER)
            w.SetMinSize((80,-1))
            GridSiz.Add(w,0,wx.ALIGN_CENTER)

        # show input for each var in expression.
        for v in self.exprVarLst:
            # label
            GridSiz.Add(wx.StaticText(self.varbox,wx.ID_ANY,v),0,wx.ALIGN_CENTER,0)
            # assignment type
            ch = wx.Choice(
                self.varbox, wx.ID_ANY,
                choices = [choices[i] for i in self.AllowedChoices]
                )
            GridSiz.Add(ch,0,wx.ALIGN_LEFT,0)
            if v in self.varSelect and self.varSelect.get(v) in self.AllowedChoices:
                i = self.AllowedChoices.index(self.varSelect[v])
                ch.SetSelection(i)
            else:
                ch.SetSelection(wx.NOT_FOUND)
            ch.label = v
            ch.Bind(wx.EVT_CHOICE,self.OnChoice)

            # var name/var assignment
            if self.varSelect.get(v) is None:
                GridSiz.Add((-1,-1),0,wx.ALIGN_LEFT|wx.EXPAND,0)
            elif self.varSelect.get(v) == 0:
                wid = G2G.ValidatedTxtCtrl(self.varbox,self.varName,v,
                                            #OnLeave=self.OnTxtLeave,
                                            size=(50,-1))
                GridSiz.Add(wid,0,wx.ALIGN_LEFT|wx.EXPAND,0)
            else:
                wid = wx.StaticText(self.varbox,wx.ID_ANY,self.varName[v])
                GridSiz.Add(wid,0,wx.ALIGN_LEFT,0)

            # value
            if self.varSelect.get(v) is None:
                GridSiz.Add((-1,-1),0,wx.ALIGN_RIGHT|wx.EXPAND,0)
            elif self.varSelect.get(v) == 0:
                wid = G2G.ValidatedTxtCtrl(self.varbox,self.varValue,v,
                                            #OnLeave=self.OnTxtLeave,
                                            size=(75,-1))
                GridSiz.Add(wid,0,wx.ALIGN_LEFT|wx.EXPAND,0)
                wid.Bind(wx.EVT_CHAR,self.OnChar)
            else:
                var = self.varName[v]
                if '*' in var:
                    #[self.parmDict[v] for v in LookupWildCard(var,self.parmDict.keys())]
                    #print self.varValue[v]
                    vs = G2obj.LookupWildCard(var,self.parmDict.keys())
                    s = '('+str(len(vs))+' values)'
                elif var in self.parmDict:
                    val = self.parmDict[var]
                    s = G2py3.FormatSigFigs(val).rstrip('0')
                else:
                    s = '?'
                wid = wx.StaticText(self.varbox,wx.ID_ANY,s)
                GridSiz.Add(wid,0,wx.ALIGN_LEFT,0)

            # show a refine flag for Free Vars only
            if self.varSelect.get(v) == 0 and self.fit:
                self.varRefflag[v] = self.varRefflag.get(v,True)
                wid = G2G.G2CheckBox(self.varbox,'',self.varRefflag,v)
                GridSiz.Add(wid,0,wx.ALIGN_LEFT|wx.EXPAND,0)
            else:
                wid = (-1,-1)
                GridSiz.Add(wid,0,wx.ALIGN_LEFT|wx.EXPAND,0)

        Siz.Add(GridSiz)
        self.varbox.SetSizer(Siz,True)
        xwid,yhgt = Siz.Fit(self.varbox)
        self.varbox.SetMinSize((xwid,130))
        self.varbox.SetAutoLayout(1)
        self.varbox.SetupScrolling()
        self.varbox.Refresh()
        self.Layout()
        #self.mainsizer.Fit(self)
        self.SendSizeEvent() # force repaint
        return

    def OnDepChoice(self,event):
        '''Respond to a selection of a variable type for a label in
        an expression
        '''
        sel = self.depChoices[event.GetEventObject().GetSelection()]
        var = self.SelectG2var(sel,'Dependent variable',self.depParmLists[sel])
        if var is None:
            self.dependentVar = None
            self.OnValidate(None)
            event.GetEventObject().SetSelection(wx.NOT_FOUND)
            return
        self.dependentVar = var
        self.depLabel.SetLabel(var)
        self.OnValidate(None)
        self.Layout()

    def GetDepVar(self):
        '''Returns the name of the dependent variable, when depVarDict is used.
        '''
        return self.dependentVar
        
    def OnChoice(self,event):
        '''Respond to a selection of a variable type for a label in
        an expression
        '''
        v = event.GetEventObject().label
        sel = self.AllowedChoices[event.GetEventObject().GetSelection()]
        if sel == 0:
            sv = G2obj.MakeUniqueLabel(v,self.usedVars)
            self.varSelect[v] = sel
            self.varName[v] = sv
            self.varValue[v] = self.varValue.get(v,0.0)
        else:
            var = self.SelectG2var(sel,v,self.parmLists[sel])
            if var is None:
                self.OnValidate(None)
                return
            self.varSelect[v] = sel
            self.varName[v] = var
        self.OnValidate(None)

    def SelectG2var(self,sel,var,parmList):
        '''Offer a selection of a GSAS-II variable. 

        :param int sel: Determines the type of variable to be selected.
          where 1 is used for Phase variables, and 2 for Histogram/Phase vars,
          3 for Histogram vars and 4 for Global vars.
        :returns: a variable name or None (if Cancel is pressed)
        '''
        if not parmList:
            return None
        l2 = l1 = 1
        for i in parmList:
            l1 = max(l1,len(i))
            loc,desc = G2obj.VarDescr(i)
            l2 = max(l2,len(loc))
        fmt = u"{:"+str(l1)+"s} {:"+str(l2)+"s} {:s}"
        varListlbl = [fmt.format(i,*G2obj.VarDescr(i)) for i in parmList]

        dlg = G2G.G2SingleChoiceDialog(
            self,'Select GSAS-II variable for '+str(var)+':',
            'Select variable',
            varListlbl,monoFont=True)
        dlg.SetSize((625,250))
        dlg.CenterOnParent()
        var = None
        if dlg.ShowModal() == wx.ID_OK:
            i = dlg.GetSelection()
            var = parmList[i]
        dlg.Destroy()
        return var

    def showError(self,msg1,msg2='',msg3=''):
        '''Show an error message of 1 to 3 sections. The second
        section is shown in an equally-spaced font. 
        
        :param str msg1: msg1 is shown in a the standard font
        :param str msg2: msg2 is shown in a equally-spaced (wx.MODERN) font
        :param str msg3: msg3 is shown in a the standard font
        '''
        self.OKbtn.Disable()
        if self.ExtraBtn: self.ExtraBtn.Disable()
        self.varSizer.Clear(True)
        self.errbox = wxscroll.ScrolledPanel(self,style=wx.HSCROLL)
        self.errbox.SetMinSize((200,130))
        self.varSizer.Add(self.errbox,1,wx.ALL|wx.EXPAND,1)
        Siz = wx.BoxSizer(wx.VERTICAL)
        errMsg1 = wx.StaticText(self.errbox, wx.ID_ANY,"")
        Siz.Add(errMsg1, 0, wx.ALIGN_LEFT|wx.LEFT|wx.EXPAND, 5)
        errMsg2 = wx.StaticText(self.errbox, wx.ID_ANY,"\n\n")
        font1 = wx.Font(errMsg2.GetFont().GetPointSize(),
                        wx.MODERN, wx.NORMAL, wx.NORMAL, False)
        errMsg2.SetFont(font1)
        Siz.Add(errMsg2, 0, wx.ALIGN_LEFT|wx.LEFT|wx.EXPAND, 5)
        errMsg3 = wx.StaticText(self.errbox, wx.ID_ANY,"")
        Siz.Add(errMsg3, 0, wx.ALIGN_LEFT|wx.LEFT|wx.EXPAND, 5)
        self.errbox.SetSizer(Siz,True)
        Siz.Fit(self.errbox)
        errMsg1.SetLabel(msg1)
        errMsg2.SetLabel("  "+msg2)
        errMsg2.Wrap(-1)
        errMsg3.SetLabel(msg3)
        self.Layout()

    def OnValidate(self,event):
        '''Respond to a press of the Validate button or when a variable
        is associated with a label (in :meth:`OnChoice`)
        '''
        self.setEvalResult('(expression cannot be evaluated)')
        self.timer.Stop()
        self.expr = self.exCtrl.GetValue().strip()
        self.varSizer.Clear(True)
        if not self.expr: 
            self.showError(
                "Invalid Expression:","",
                "(an expression must be entered)")
            return
        exprObj = G2obj.ExpressionObj()
        ret = exprObj.ParseExpression(self.expr)
        if not ret:
            self.showError(*exprObj.lastError)
            return
        self.exprVarLst,pkgdict = ret
        wx.CallAfter(self.Repaint,exprObj)
            
    def Repaint(self,exprObj):
        'Redisplay the variables and continue the validation'
        self.ShowVars() # show widgets to set vars
        msg = self.CheckVars() 
        if msg:
            self.setEvalResult(msg)
            return
        exprObj.LoadExpression(
            self.expr,
            self.exprVarLst,
            self.varSelect,
            self.varName,
            self.varValue,
            self.varRefflag,
            )
        try:
            calcobj = G2obj.ExpressionCalcObj(exprObj)
            calcobj.SetupCalc(self.parmDict)
            val = calcobj.EvalExpression()
        except Exception as msg:
            self.setEvalResult("Error in evaluation: "+str(msg))
            return
        if not np.isfinite(val):
            self.setEvalResult("Expression value is infinite or out-of-bounds")
            return
        s = G2py3.FormatSigFigs(val).rstrip('0')
        depVal = ""
        if self.depVarDict:
            if not self.dependentVar:
                self.setEvalResult("A dependent variable must be selected.")
                return
            depVal = '; Variable "' + self.dependentVar + '" = ' + str(
                self.depVarDict.get(self.dependentVar,'?')
                )
        self.setEvalResult("Expression evaluates to: "+str(s)+depVal)
        self.OKbtn.Enable()
        if self.ExtraBtn: self.ExtraBtn.Enable()
        
if __name__ == "__main__":
    app = wx.PySimpleApp() # create the App
    frm = wx.Frame(None)
    frm.Show()
    PSvarDict = {'::a':1.0,'::b':1.1,'0::c':1.2,'::AlongVaraiableName':1000.}
    #PSvars = PSvarDict.keys()
    indepvarDict = {'Temperature':1.0,'Pressure':1.1,'Phase of Moon':1.2,'1:1:HAP':1.3}
    dlg = ExpressionDialog(frm,indepvarDict,
                           header="Edit the PseudoVar expression",
                           fit=False,
                           depVarDict=PSvarDict,
                           #VarLabel="New PseudoVar",                           
                           )
    newobj = dlg.Show(True)
    dlg = ExpressionDialog(frm,indepvarDict,newobj,
                           header="Edit the PseudoVar expression",
                           fit=False,
                           depVarDict=PSvarDict,
                           #VarLabel="New PseudoVar",                           
                           )
    newobj = dlg.Show(True)
    print dlg.GetDepVar()
    dlg = ExpressionDialog(frm,PSvarDict,
                           header="Edit the PseudoVar expression",
                           fit=True)
    newobj = dlg.Show(True)
    print dlg.GetDepVar()
    import sys
    #sys.exit()

    #app.MainLoop()


    import cPickle
    def showEQ(calcobj):
        print
        print calcobj.eObj.expression
        for v in sorted(calcobj.eObj.freeVars.keys()+calcobj.eObj.assgnVars.keys()):
            print "  ",v,'=',calcobj.exprDict[v]
        print calcobj.EvalExpression()
    print "starting test"
    obj = G2obj.ExpressionObj()
    obj.expression = "A*np.exp(B)"
    obj.assgnVars =  {'B': '0::Afrac:*'}
    obj.freeVars =  {'A': [u'A', 0.5, True]}
    obj.CheckVars()
    parmDict2 = {'0::Afrac:0':1.0, '0::Afrac:1': 1.0}
    calcobj = G2obj.ExpressionCalcObj(obj)
    calcobj.SetupCalc(parmDict2)
    showEQ(calcobj)
    fp = open('/tmp/obj.pickle','w')
    cPickle.dump(obj,fp)
    fp.close()
    
    obj.expression = "A*np.exp(-2/B)"
    obj.assgnVars =  {'A': '0::Afrac:0', 'B': '0::Afrac:1'}
    obj.freeVars =  {}
    parmDict1 = {'0::Afrac:0':1.0, '0::Afrac:1': -2.0}
    calcobj = G2obj.ExpressionCalcObj(obj)
    calcobj.SetupCalc(parmDict1)
    showEQ(calcobj)

    fp = open('/tmp/obj.pickle','r')
    obj = cPickle.load(fp)
    fp.close()
    parmDict2 = {'0::Afrac:0':0.0, '0::Afrac:1': 1.0}
    calcobj = G2obj.ExpressionCalcObj(obj)
    calcobj.SetupCalc(parmDict2)
    showEQ(calcobj)

    parmDict2 = {'0::Afrac:0':1.0, '0::Afrac:1': 1.0}
    calcobj = G2obj.ExpressionCalcObj(obj)
    calcobj.SetupCalc(parmDict2)
    showEQ(calcobj)
    
