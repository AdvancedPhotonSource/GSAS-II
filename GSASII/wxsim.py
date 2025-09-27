'''A library of routines used to simulate user actions for 
playback of GSAS-II actions. 
'''

import time
import wx
from GSASII.GSASIIdataGUI import GetGPXtreeItemId as GetGPXtreeItemId

asyncResults = {}
# generic interface routines to get wx objects from widget names
# and interact with those widgets
def getAllChildren(parent):
    '''Recursively finds all children of a given parent window.

    :Returns: A list containing all descendant wx.Window objects.
    '''
    all_children = []
    for child in parent.GetChildren():
        all_children.append(child)
        all_children.extend(getAllChildren(child))
    return all_children

def getButtons(winname, btnList):
    '''Finds a window with a label matching `winname` and
    then the button objects matching the names in `btnlist`

    If run in a thread, this should be called by CallAfter
    '''    
    asyncResults.clear()
    # find window
    for w in wx.GetTopLevelWindows():
        if w.GetTitle() == winname:
            #print(f"Window {winname!r} found")
            win = w
            asyncResults['frame'] = win
            break
    else:
        asyncResults['error'] = f"Window {winname!r} not found"        
        return
    # find each button
    btns = []
    for btnname in btnList:
        for b in getAllChildren(win):
            if not hasattr(b,'GetLabelText'): continue # not button
            if not issubclass(b.__class__,wx.Button): continue # not button
            if b.GetLabelText() == btnname:
                btns.append(b)
                break
        else:
            asyncResults['error'] = f'Button {btnname} not found'
            return
    asyncResults['buttons'] = btns

def invokeButton(frame, btn):
    '''"press" button (`btn`) in a supplied window (`frame`) by 
    moving the mouse to the center of the button and "left-clicking"
    via UIActionSimulator.

    If run in a thread, this should be called by CallAfter
    '''
    pos = frame.ClientToScreen(btn.GetPosition() + btn.GetSize()/2)
    sim = wx.UIActionSimulator()
    sim.MouseMove(pos.x,pos.y)
    time.sleep(0.05)  # this should be short so that there is
    # no chance for the user to move the mouse, but long enough so that the
    # move registers
    sim.MouseClick(wx.MOUSE_BTN_LEFT)
    #print('invoke done @',pos)
    time.sleep(0.2)
    asyncResults['done'] = True

def pressButtonsInWindow(winname,bnlList):
    '''Run this in a thread to "press buttons" in a dialog
    that is opened just after this is started. Note that 
    final button should close the dialog.

    This will need more development to do more things in a
    "semi-modal" dialog
    '''
    time.sleep(0.25) # wait for the window to open
    
    asyncResults.clear()
    wx.CallAfter(getButtons,winname,bnlList)
    count = waitForDone(True)
    print('getButtons',asyncResults, count)
    buttons,frame = asyncResults['buttons'],asyncResults['frame']
    #wx.CallAfter(wx.Yield)
    
    for i, b in enumerate(buttons):
        asyncResults.clear()
        wx.CallAfter(invokeButton,frame,b)
        print(f'button {i} pressed')
        count = waitForDone(True)
        #if i < len(buttons) - 1:  # If not the last button
        #    time.sleep(0.3)  # Wait between button clicks
        time.sleep(0.3)  # Wait between button clicks
    
def GetPhaseItem(G2frame,findname=None,findcount=None):
    '''Select a phase from the data tree by name or by count.
    Specify a name for the phase or the number of the phase.

    :param wx.Frame G2frame: reference to main GSASII frame
    :param str findname: name of phase to be found (or None to 
       use findcount)
    :param str findcount: index of phase to be found (or None to 
       use findname). Starts with 0.

    This is intended to be called from wx.CallAfter(),
    so results are placed in dict asyncResults which
    must be defined before this is called:
     * asyncResults['count'] has the index for the located phase 
     * asyncResults['name']  has the name for the located phase
     * asyncResults['error'] has error information, if an error occurs
     * asyncResults['done'] if findcount >= number of phase entries

    '''
    asyncResults.clear()
    if findname is None and findcount is None:
        #print('findname=None & findcount=None')
        asyncResults['error'] = 'invalid GetPhaseItem call'
        return
    phId = GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
    if not phId:
        #print('no phases')
        asyncResults['error'] = 'no phases'
        return
    # cycle through phases
    item, cookie = G2frame.GPXtree.GetFirstChild(phId)
    count = 0
    while item:
        name = G2frame.GPXtree.GetItemText(item)
        #print('Phase',name,count)
        if findname is not None and findname == name:
            break
        if findcount is not None and findcount == count:
            break
        item, cookie = G2frame.GPXtree.GetNextChild(phId, cookie)
        count += 1
    else:
        if findname is not None:
            asyncResults['error'] = f'phase {findname!r} not found'
        else:
            asyncResults['done'] = f'phase #{findcount} not found'
        return
    G2frame.GPXtree.SelectItem(item)
    asyncResults['count'] = count
    asyncResults['name'] = name

def GetTab(G2frame,findname):
    '''Select a tab on the current data window pane

    :param wx.Frame G2frame: reference to main GSASII frame
    :param str findname: name of tab to be found

    This is intended to be called from wx.CallAfter(),
    so results are placed in dict asyncResults which
    must be defined before this is called:
     * asyncResults['error'] has error information, if an error occurs

    '''
    asyncResults.clear()
    for i in range(G2frame.phaseDisplay.GetPageCount()):
        if G2frame.phaseDisplay.GetPageText(i) == findname:
            G2frame.phaseDisplay.SetSelection(i)
            asyncResults['done'] = True
            return
    else:
        asyncResults['error'] = f'tab {finaname!r} not found'
        #print(asyncResults['error'])
        return

def InvokeMenuCommand(G2frame,menuname,itemname):
    '''Finds and invokes a menu command from the currently displayed menu

    :param wx.Frame G2frame: reference to main GSASII frame
    :param str menuname: name of menu (File, Help,...)
    :param str itemname: name of menu (Open, Quit,...)

    This is intended to be called from wx.CallAfter(),
    so results are placed in dict asyncResults which
    must be defined before this is called:
     * asyncResults['menuitem'] has the menu object for the selected 
     * asyncResults['error'] has error information, if an error occurs

    '''
    asyncResults.clear()
    # select menu command
    for i in range(G2frame.GetMenuBar().GetMenuCount()):
        if G2frame.GetMenuBar().GetMenuLabelText(i) == menuname:
            menu = G2frame.GetMenuBar().GetMenu(i)
            break
    else:
        asyncResults['error'] = f'menu {menuname!r} not found'
        return
    for i in menu.GetMenuItems():
        if menu.GetLabelText(i.Id) == itemname:
            menuitem = i
            break
    else:
        asyncResults['error'] = f'menu command {itemname!r} not found'
        return
    
    try:
        menuevent = wx.CommandEvent(wx.EVT_MENU.typeId, menuitem.Id)
        G2frame.ProcessEvent(menuevent)
    except Exception as msg:
        asyncResults['error'] = f'Error from InvokeMenuCommand: {msg}'
    asyncResults['menuitem'] = menuitem

def SelectRowsInPhaseGrid(G2frame,selections='all'):
    '''This selects rows in a grid on the currently displayed 
    phase notebook
    '''
    # find grid -- assume only one
    for i in getAllChildren(G2frame.phaseDisplay.GetCurrentPage()):
        if i.__class__.__base__ is wx._grid.Grid:
            grid = i
            break
    else:
        asyncResults['error'] = 'Error no grid found'
        return
    if selections == 'all':
        selections = list(range(grid.GetNumberRows()))
    for i in range(grid.GetNumberRows()):
        if i in selections: grid.SelectRow(i,True)
    asyncResults['done'] = True
            
def waitForDone(verbose=False):
    time.sleep(0.1)
    count = 0
    while len(asyncResults) == 0:
        count += 1
        time.sleep(0.05)
        if count > 1000: raise Exception('too many waits')
    if 'error' in asyncResults:
        raise Exception(f"async failure: {asyncResults['error']}")
    if verbose: print('done',count,asyncResults)
    return count

# def waitForModal(verbose=False):
#     time.sleep(0.1)
#     count = 0
#     while len(asyncResults) == 0:
#         count += 1
#         time.sleep(0.05)
#         if count > 1000: raise Exception('too many waits')
#     if 'error' in asyncResults:
#         raise Exception(f"async failure: {asyncResults['error']}")
#     if verbose: print('done',count,asyncResults)
#     return count
