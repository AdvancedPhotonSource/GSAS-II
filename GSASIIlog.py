# -*- coding: utf-8 -*-
#GSASIIlog - Routines used to track and replay "actions"
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*GSASIIlog: Logging of "Actions"*
---------------------------------

Module to provide logging services, e.g. track and replay "actions"
such as menu item, tree item, button press, value change and so on. 
'''
from __future__ import division, print_function
try:
    import wx
except ImportError:
    pass
import GSASIIdataGUI as G2gd
import GSASIIpath
# Global variables 
MenuBindingLookup = {}
'Lookup table for Menu buttons'
ButtonBindingLookup = {}
'Lookup table for button objects'
G2logList = [None]
'Contains a list of logged actions; first item is ignored'
LogInfo = {'Logging':False, 'Tree':None, 'LastPaintAction':None}
'Contains values that are needed in the module for past actions & object location'

# TODO:
# Might want to save the last displayed screen with some objects to make sure
# the commands are executed in a sensible order
# 1) save tree press under tab press item.
# 2) save tree & tab for button press

# TODO:
# probably need an object for saving calls and arguments to call specific functions/methods.
# The items to be called need to register themselves so that they can be found
#   This needs to be done with a Register(function,'unique-string') call after every def
#   and then a LogCall('unique-string',pargs,kwargs) call inside the routine

debug = GSASIIpath.GetConfigValue('logging_debug')
#===========================================================================
# objects for logging variables
def _l2s(lst,separator='+'):
    'Combine a list of objects into a string, with a separator'
    s = ''
    for i in lst: 
        if s != '': s += separator
        s += '"'+str(i)+'"'
    return s

class LogEntry(object):
    '''Base class to define logging objects. These store information on events
    in a manner that can be pickled and saved -- direct references to wx objects
    is not allowed.

    Each object must define:
    
     *  __init__: stores the information needed to log & later recreate the action 
     *  __str__ : shows a nice ASCII string for each action
     *  Replay:   recreates the action when the log is played
     
    optional:
    
      * Repaint:  redisplays the current window
      
    '''
    def __init__(self):
        # Must be overridden
        raise Exception('No __init__ defined')
    def __str__(self):
        # Must be overridden
        raise Exception('No __str__ defined')
    def Replay(self):
        # Must be overridden
        raise Exception('No Replay defined')
    def Repaint(self):
        # optional
        pass

class VarLogEntry(LogEntry):
    'object that tracks changes to a variable'
    def __init__(self,treeRefs,indexRefs,value):
        self.treeRefs = treeRefs
        self.indexRefs = indexRefs
        self.value = value
        if debug: print ('Logging var change: w/treeRefs',treeRefs,'indexRefs',indexRefs,'new value=',value)
    def __str__(self):
        treeList = self.treeRefs[:]
        if type(treeList[0]) is tuple:
            treeList[0] = 'Hist #'+str(treeList[0][1]+1)+' of type '+treeList[0][0]
        elif len(treeList) > 1 and type(treeList[1]) is int:
            treeList[1] = 'Phase #'+str(treeList[1]+1)
        return 'Variable change: Key(s)= '+_l2s(self.indexRefs)+' to value='+str(self.value)
    def Replay(self):
        'Perform a Variable Change action, when read from the log'
        parentId = LogInfo['Tree'].root
        for i,treeitem in enumerate(self.treeRefs):
            if i == 0 and type(treeitem) is tuple:
                treeitem = LogInfo['Tree'].ConvertRelativeHistNum(*treeitem)
            item, cookie = LogInfo['Tree'].GetFirstChild(parentId)
            while item:
                if LogInfo['Tree'].GetItemText(item) == treeitem:
                    parentId = item
                    break
                else:
                    item, cookie = LogInfo['Tree'].GetNextChild(parentId, cookie)
            else:
                raise Exception("Tree item not found for "+str(self))
        # get the inner most data array
        data = LogInfo['Tree'].GetItemPyData(item)
        for item in self.indexRefs[:-1]:
            data = data[item]
        # set the value
        data[self.indexRefs[-1]] = self.value

class MenuLogEntry(LogEntry):
    'object that tracks when a menu command is executed'
    def __init__(self,menulabellist):
        self.menulabellist = menulabellist
        if debug:
            t = menulabellist[:]
            t.reverse()
            l = ''
            for item in t:
                if l: l += ' -> '
                l += item
            print ('Logging menu command: '+l)
    def __str__(self):
        return 'Menu press: From '+_l2s(self.menulabellist,'/')
    def Replay(self):
        'Perform a Menu item action when read from the log'
        key = ''
        for item in self.menulabellist:
            if key: key += '+'
            key += item
        if MenuBindingLookup.get(key):
            handler,id,menuobj = MenuBindingLookup[key]
            MyEvent = wx.CommandEvent(wx.EVT_MENU.typeId, id)
            MyEvent.SetEventObject(menuobj)
            handler(MyEvent)
        else:
            raise Exception('No binding for menu item '+key)        
            
class TabLogEntry(LogEntry):
    'Object to track when tabs are pressed in the DataFrame window'
    def __init__(self,title,tabname):
        self.wintitle = title
        self.tablabel = tabname
        if debug: print ('Logging tab: "'+tabname+'" on window titled '+title)
    def __str__(self):
        return 'Tab press: Tab='+_l2s([self.tablabel])+' on window labeled '+str(self.wintitle)
    def Repaint(self):
        'Used to redraw a window created in response to a Tab press'
        if debug: print ('Repaint')
        saveval = LogInfo['LastPaintAction']
        self.Replay()
        LogInfo['LastPaintAction'] = saveval
    def Replay(self):
        'Perform a Tab press action when read from the log'
        wintitle = self.wintitle
        tabname = self.tablabel
        LogInfo['LastPaintAction'] = self
        if LogInfo['Tree'].G2frame.GetTitle() != wintitle:
            print (LogInfo['Tree'].G2frame.GetTitle()+' != '+wintitle)
            raise Exception('tab in wrong window')
        for PageNum in range(LogInfo['Tree'].G2frame.dataDisplay.GetPageCount()):
            if tabname == LogInfo['Tree'].G2frame.dataDisplay.GetPageText(PageNum):
                LogInfo['Tree'].G2frame.dataDisplay.SetSelection(PageNum)
                return
        else:
            print (tabname+'not in'+[
                LogInfo['Tree'].G2frame.dataDisplay.GetPageText(PageNum) for
                PageNum in range(LogInfo['Tree'].G2frame.dataDisplay.GetPageCount())])
            raise Exception('tab not found')
def MakeTabLog(title,tabname):
    'Create a TabLogEntry action log'
    G2logList.append(TabLogEntry(title,tabname))

class TreeLogEntry(LogEntry):
    'Object to track when tree items are pressed in the main window'
    def __init__(self,itemlist):
        self.treeItemList = itemlist
        if debug: print ('Logging press on tree: '+itemlist)
    def __str__(self):
        treeList = self.treeItemList[:]
        if type(treeList[0]) is tuple:
            treeList[0] = 'Hist #'+str(treeList[0][1]+1)+' of type '+treeList[0][0]
        elif len(treeList) > 1 and type(treeList[1]) is int:
            treeList[1] = 'Phase #'+str(treeList[1]+1)
        return 'Tree item pressed: '+_l2s(treeList)
    def Repaint(self):
        'Used to redraw a window created in response to a click on a data tree item'
        if debug: print ('Repaint')
        saveval = LogInfo['LastPaintAction']
        LogInfo['Tree'].SelectItem(LogInfo['Tree'].root) # need to select something else
        wx.Yield()
        self.Replay()
        LogInfo['LastPaintAction'] = saveval
    def Replay(self):
        'Perform a Tree press action when read from the log'
        LogInfo['LastPaintAction'] = self
        parent = LogInfo['Tree'].root
        for i,txt in enumerate(self.treeItemList):
            if i == 0 and type(txt) is tuple:
                txt = LogInfo['Tree'].ConvertRelativeHistNum(*txt)
            elif i == 1 and type(txt) is int and self.treeItemList[0] == "Phases":
                txt = LogInfo['Tree'].ConvertRelativePhaseNum(txt)
            item = G2gd.GetGPXtreeItemId(LogInfo['Tree'].G2frame,parent,txt)
            if not item:
                print ('Not found',+txt)
                return
            else:
                parent = item
        else:
            LogInfo['Tree'].SelectItem(item)
def MakeTreeLog(textlist):
    'Create a TreeLogEntry action log'
    G2logList.append(TreeLogEntry(textlist))
    
class ButtonLogEntry(LogEntry):
    'Object to track button press'
    def __init__(self,locationcode,label):
        self.locationcode = locationcode
        self.label = label
        if debug: print ('Logging '+label+' button press in '+locationcode)
    def __str__(self):
        return 'Press of '+self.label+' button in '+self.locationcode
    def Replay(self):
        key = self.locationcode + '+' + self.label
        if ButtonBindingLookup.get(key):
            btn = ButtonBindingLookup[key]
            clickEvent = wx.CommandEvent(wx.EVT_BUTTON.typeId, btn.GetId())
            clickEvent.SetEventObject(btn)
            #btn.GetEventHandler().ProcessEvent(clickEvent)
            btn.handler(clickEvent)
def MakeButtonLog(locationcode,label):
    'Create a ButtonLogEntry action log'
    G2logList.append(ButtonLogEntry(locationcode,label))
    
            
def _wrapper(func):
    def _wrapped(self, *args, **kwargs):
        return getattr(self.obj, func)(*args, **kwargs)
    return _wrapped

class DictMeta(type):
    def __new__(cls, name, bases, dct):
        default_attrs = dir(object) + ['__getitem__', '__str__']
        for attr in dir(dict):
            if attr not in default_attrs:
                dct[attr] = _wrapper(attr)
        return type.__new__(cls, name, bases, dct)

class dictLogged(object):
    '''A version of a dict object that tracks the source of the
    object back to the location on the G2 tree.
    If a list (tuple) or dict are pulled from inside this object
    the source information is appended to the provinance tracking
    lists.

    tuples are converted to lists.
    '''
    __metaclass__ = DictMeta

    def __init__(self, obj, treeRefs, indexRefs=[]):
        self.treeRefs = treeRefs
        self.indexRefs = indexRefs
        self.obj = obj

    def __getitem__(self,key):
        val = self.obj.__getitem__(key)   
        if type(val) is tuple:
            #if debug: print 'Converting to list',key
            val = list(val)
            self.obj[key] = val
        if type(val) is dict:
            #print 'dict'
            return dictLogged(val,self.treeRefs,self.indexRefs+[key])
        elif type(val) is list:
            #print 'list'
            return listLogged(val,self.treeRefs,self.indexRefs+[key])
        else:
            #print type(val)
            return val

    def __str__(self):
        return self.obj.__str__()
    def __repr__(self):
        return self.obj.__str__() + " : " + str(self.treeRefs) + ',' + str(self.indexRefs)

class ListMeta(type):
    def __new__(cls, name, bases, dct):
        default_attrs = dir(object) + ['__getitem__', '__str__']
        for attr in dir(list):
            if attr not in default_attrs:
                dct[attr] = _wrapper(attr)
        return type.__new__(cls, name, bases, dct)

class listLogged(object):
    '''A version of a list object that tracks the source of the
    object back to the location on the G2 tree.
    If a list (tuple) or dict are pulled from inside this object
    the source information is appended to the provinance tracking
    lists.
    
    tuples are converted to lists.
    '''
    __metaclass__ = ListMeta

    def __init__(self, obj, treeRefs, indexRefs=[]):
        self.treeRefs = treeRefs
        self.indexRefs = indexRefs
        self.obj = obj

    def __getitem__(self,key):
        val = self.obj.__getitem__(key)   
        if type(val) is tuple:
            #if debug: print 'Converting to list',key
            val = list(val)
            self.obj[key] = val
        if type(val) is dict:
            #print 'dict'
            return dictLogged(val,self.treeRefs,self.indexRefs+[key])
        elif type(val) is list:
            #print 'list'
            return listLogged(val,self.treeRefs,self.indexRefs+[key])
        else:
            #print type(val)
            return val

    def __str__(self):
        return self.obj.__str__()
    def __repr__(self):
        return self.obj.__str__() + " : " + str(self.treeRefs) + ',' + str(self.indexRefs)

#===========================================================================
# variable tracking
def LogVarChange(result,key):
    'Called when a variable is changed to log that action'
    if not LogInfo['Logging']: return
    if hasattr(result,'treeRefs'):
        treevars = result.treeRefs[:]
        treevars[0] = LogInfo['Tree'].GetRelativeHistNum(treevars[0])
        if treevars[0] == "Phases" and len(treevars) > 1:
            treevars[1] = LogInfo['Tree'].GetRelativePhaseNum(treevars[1])
        lastLog = G2logList[-1]
        fullrefs = result.indexRefs+[key]
        if type(lastLog) is VarLogEntry:
            if lastLog.treeRefs == treevars and lastLog.indexRefs == fullrefs:
                lastLog.value = result[key]
                if debug: print ('update last log to '+result[key])
                return
        G2logList.append(VarLogEntry(treevars,fullrefs,result[key]))
    else:
        print (key+' Error: var change has no provenance info')

#===========================================================================
# menu command tracking
def _getmenuinfo(id,G2frame,handler):
    '''Look up the menu/menu-item label tree from a menuitem's Id
    
    Note that menubars contain multiple menus which contain multiple menuitems.
    A menuitem can itself point to a menu and if so that menu can contain 
    multiple menuitems.

    Here we start with the last menuitem and look up the label for that as well
    as all parents, which will be found in parent menuitems (if any) and the menubar.
    
        menuitem    ->  menu ->  menubar
           |                          |
           |->itemlabel               |-> menulabel
           
    or

        menuitem    ->  (submenu -> menuitem)*(n times) -> menu -> menubar
           |                            |                             |
           |->itemlabel                 |-> sublabel(s)               |-> menulabel
           
    :returns: a list containing all the labels and the menuitem object
       or None if the menu object will not be cataloged. 
    '''
    # don't worry about help menuitems
    if id == wx.ID_ABOUT: return
    # get the menu item object by searching through all menubars and then its label
    for menubar in G2frame.dataMenuBars:
        menuitem = menubar.FindItemById(id)
        if menuitem:
            #print 'getmenuinfo found',id,menuitem
            break
    else:
        #print '****** getmenuinfo failed for id=',id,'binding to=',handler
        #raise Exception('debug: getmenuinfo failed')
        return
    menuLabelList = [menuitem.GetItemLabel()]
    
    # get the menu where the current item is located
    menu = menuitem.GetMenu()
    while menu.GetParent(): # is this menu a submenu of a previous menu?
        parentmenu = menu.GetParent()
        # cycle through the parentmenu until we find the menu
        for i in range(parentmenu.GetMenuItemCount()):
            if parentmenu.FindItemByPosition(i).GetSubMenu()==menu:
                menuLabelList += [parentmenu.FindItemByPosition(i).GetItemLabel()]
                break
        else:
            # menu not found in menu, something is wrong
            print ('error tracing menuitem to parent menu'+menuLabelList)
            #raise Exception('debug1: error tracing menuitem')
            return
        menu = parentmenu
        
    i,j= wx.__version__.split('.')[0:2]
    if int(i)+int(j)/10. > 2.8 and 'wxOSX' in wx.PlatformInfo:
        # on mac, with wx 2.9+ the menubar has a menu and this is found above, so
        # we are now done.
        return menuLabelList,menuitem
    
    menubar = menu.MenuBar
    for i in range(menubar.GetMenuCount()):
        if menubar.GetMenu(i) == menu:
            menuLabelList += [menubar.GetMenuLabel(i)]
            return menuLabelList,menuitem

    # menu not found in menubar, something is wrong
    print ('error tracing menuitem to menubar'+menuLabelList)
    #raise Exception('debug2: error tracing menuitem')
    return

def SaveMenuCommand(id,G2frame,handler):
    '''Creates a table of menu items and their pseudo-bindings
    '''
    menuinfo = _getmenuinfo(id,G2frame,handler)
    if not menuinfo: return
    menuLabelList,menuobj = menuinfo
    key = ''
    for item in menuLabelList:
        if key: key += '+'
        key += item
    MenuBindingLookup[key] = [handler,id,menuobj]
    return menuLabelList

def InvokeMenuCommand(id,G2frame,event):
    '''Called when a menu item is used to log the action as well as call the
    routine "bind"ed to that menu item
    '''
    menuLabelList,menuobj = _getmenuinfo(id,G2frame,None)
    key = ''
    if menuLabelList: 
        for item in menuLabelList:
            if key: key += '+'
            key += item
    if key in MenuBindingLookup:
        if LogInfo['Logging']: 
            G2logList.append(MenuLogEntry(menuLabelList))
        handler = MenuBindingLookup[key][0]
        handler(event)
    else:
        print ('Error no binding for menu command'+menuLabelList+'id=%d'%id)
        return

#===========================================================================
# Misc externally callable routines
def LogOn():
    'Turn On logging of actions'
    if debug: print ('LogOn')
    LogInfo['Logging'] = True

def LogOff():
    'Turn Off logging of actions'
    if debug: print ('LogOff')
    LogInfo['Logging'] = False
    
def ShowLogStatus():
    'Return the logging status'
    return LogInfo['Logging']

def OnReplayPress(event):
    'execute one or more commands when the replay button is pressed'
    clb = LogInfo['clb']
    dlg = clb.GetTopLevelParent()
    sels = sorted(clb.GetSelections())
    if not sels:
        dlg1 = wx.MessageDialog(dlg,
            'Select one or more items in the list box to replay',
            'No selection actions',
            wx.OK)
        dlg1.CenterOnParent()
        dlg1.ShowModal()
        dlg1.Destroy()
        return
    logstat = ShowLogStatus()
    if logstat: LogOff()
    if debug: print (70*'=')
    for i in sels:
        i += 1
        item = G2logList[i]
        if debug: print ('replaying'+item)
        item.Replay()
        wx.Yield()
    if i >= len(G2logList)-1:
        dlg.EndModal(wx.ID_OK)
    else:
        clb.DeselectAll()
        clb.SetSelection(i)
    if debug: print (70*'=')
    # if the last command did not display a window, repaint it in
    # case something on that window changed. 
    if item != LogInfo['LastPaintAction'] and hasattr(LogInfo['LastPaintAction'],'Repaint'):
        LogInfo['LastPaintAction'].Repaint()
    if logstat: LogOn()
     
def ReplayLog(event):
    'replay the logged actions'
    LogInfo['LastPaintAction'] = None # clear the pointed to the last data window
    # is this really needed? -- probably not. 
    commandList = []
    for item in G2logList:
        if item: # skip over 1st item in list (None)
            commandList.append(str(item))
    if not commandList:
        dlg = wx.MessageDialog(LogInfo['Tree'],
            'No actions found in log to replay',
            'Empty Log',
            wx.OK)
        dlg.CenterOnParent()
        dlg.ShowModal()
        dlg.Destroy()
        return
    dlg = wx.Dialog(LogInfo['Tree'],wx.ID_ANY,'Replay actions from log',
        style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.CENTRE)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5))
    clb = wx.ListBox(dlg, wx.ID_ANY, (30,100), wx.DefaultSize, commandList,
                     style=wx.LB_EXTENDED)
    LogInfo['clb'] = clb
    mainSizer.Add(clb,1,wx.EXPAND,1)
    mainSizer.Add((5,5))
    btn = wx.Button(dlg, wx.ID_ANY,'Replay selected')
    btn.Bind(wx.EVT_BUTTON,OnReplayPress)
    mainSizer.Add(btn,0,wx.ALIGN_CENTER,0)
    btnsizer = wx.StdDialogButtonSizer()
    OKbtn = wx.Button(dlg, wx.ID_OK,'Close')
    #OKbtn = wx.Button(dlg, wx.ID_CLOSE)
    OKbtn.SetDefault()
    OKbtn.Bind(wx.EVT_BUTTON,lambda event: dlg.EndModal(wx.ID_OK))
    btnsizer.AddButton(OKbtn)
    btnsizer.Realize()
    mainSizer.Add((-1,5),1,wx.EXPAND,1)
    mainSizer.Add(btnsizer,0,wx.ALIGN_CENTER,0)
    mainSizer.Add((-1,5))
    dlg.SetSizer(mainSizer)
    dlg.CenterOnParent()
    clb.SetSelection(0)
    dlg.ShowModal()
    dlg.Destroy()
    LogInfo['Tree'].G2frame.OnMacroRecordStatus(None) # sync the menu checkmark(s)
    return

if debug: LogOn() # for debug, start with logging enabled
