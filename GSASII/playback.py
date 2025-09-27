import time
import threading
import wx
import wxsim
from importlib import reload
reload(wxsim)
#from GSASII.GSASIIdataGUI import GetGPXtreeItemId as GetGPXtreeItemId

pause = lambda delay: time.sleep(delay)

asyncResults = {}
def Playback(G2frame):
    '''This is a "demo" routine that loads the first and then the 
    second phase from a .gpx file

    This will be run from a command callback so GUI routines can be called
    directly.
    '''
    i = 0
    G2frame.Playback = True
    try:
        for count in range(2): # cycle through two phases
            wxsim.GetPhaseItem(G2frame,None,count)
            wx.Yield()
            wxsim.waitForDone()
            wxsim.GetTab(G2frame,'Draw Atoms')
            wx.Yield()
            wxsim.waitForDone()
            threading.Thread(target=wxsim.pressButtonsInWindow,
                        args=('Select atoms for action',['Set All','OK']),
                        ).start()
            wxsim.InvokeMenuCommand(G2frame,'Edit Figure','Fill unit cell')
            wx.Yield()
            #wxsim.waitForDone(True)
            #pause(2.1)
            #print(wxsim.asyncResults)
            #pause(0.1)
        # split the window
        #G2frame.G2plotNB.nb.Split(0,wx.Right)
    finally:
    # done, reset the playback flag
        G2frame.Playback = False




# asyncResults = {}
# def Playback(G2frame):
#     '''This is a "demo" routine that loads the first and then the 
#     second phase from a .gpx file

#     This will be run in a thread, so it should use CallAfter for 
#     any actions that call wx routines.
#     '''
#     G2frame.Playback = True
#     try:
#         for count in range(2): # cycle through two phases
#             wx.CallAfter(wxsim.GetPhaseItem,G2frame,None,count)
#             wxsim.waitForDone()
#             wx.CallAfter(wxsim.GetTab,G2frame,'Draw Atoms')
#             wxsim.waitForDone()
#             wx.CallAfter(wxsim.GetTab,G2frame,'Draw Atoms')
#             wxsim.waitForDone()
#             # select all atoms
#             wx.CallAfter(wxsim.SelectRowsInPhaseGrid,G2frame)
#             # fill the cell
#             wxsim.waitForDone()
#             wx.CallAfter(wxsim.InvokeMenuCommand,G2frame,'Edit Figure','Fill unit cell')
#             wxsim.waitForDone()
#             #pause(2.1)
#             #print(wxsim.asyncResults)
#             #pause(0.1)
#             # split the window
#         wx.CallAfter(G2frame.G2plotNB.nb.Split,0,wx.Right)
#     finally:
#     # done, reset the playback flag
#         G2frame.Playback = False
#     #import GSASII.GSASIIpath
#     #GSASII.GSASIIpath.IPyBreak()
    
# def doAfterFileLoad(G2frame): # not exactly logging playback, but a start
#     # runs 
#     pause(0.2)
#     # get Phase tab
#     launch = wx.CallAfter
#     #def launch(function,*args,**kwargs):
#     #    function(*args,**kwargs)
#     #    wx.Yield()
#     for count in range(2): # cycle through two phases
#         launch(GetPhaseItem,G2frame,None,count)
#         pause(0.5)
#         if 'error' in asyncResults: return False
#         if 'done' in asyncResults: break
#         #print(asyncResults)
#         launch(GetTab,G2frame,'Draw Atoms')
#         pause(0.1)
#         if 'error' in asyncResults: return False
#         #print(asyncResults)
#         launch(GetMenuCommand,G2frame,'Edit Figure','Fill unit cell')
#         pause(0.1)
#         if 'error' in asyncResults: return False
#         #print(asyncResults)
#         # menuitem = 
#         # launch(InvokeMenuCommand,G2frame,menuitem,
#         #                  pressButtonsInWindow,
#         #                  ('Select atoms for action',('Set All','OK')))

#         # launch thread to run commands in created window
#         threading.Thread(target=doInWindow,
#             args=('Select atoms for action',['Set All','OK']),
#             #kwargs=None
#             ).start()
#         # open window. command above should close it
#         launch(InvokeMenuCommand,G2frame,asyncResults['menuitem'])
#         return
#         #breakpoint()



# def InvokeMenuCommand(G2frame,menuitem,concurrentCommand=None,
#                           concurrentArgs=None, concurrentKWargs=None):
#     '''Executes a menu found in :func:`GetMenuCommand`, optionally 
#     running a set of commands in a separate thread (needed if
#     `menuitem` opens a modal dialog.

#     :param wx.Frame G2frame: reference to main GSASII frame
#     :param wx.MenuItem menuitem: reference to a menu item
#     :param function concurrentCommand: if specified, a callable
#        routine that will be run in a separate thread while `menuitem`
#        is run.
#     :param list concurrentArgs: if specified, a list of positional 
#        argument(s) to be passed to routine `concurrentCommand`
#     :param dict concurrentKWargs: if specified, a dict of keyword 
#        argument(s) to be passed to routine `concurrentCommand`

#     This is intended to be called from wx.CallAfter(),
#     so results are placed in dict asyncResults which
#     must be defined before this is called:
#      * asyncResults['error'] has error information, if an error occurs

#     '''
#     asyncResults.clear()
#     try:
#         menuevent = wx.CommandEvent(wx.EVT_MENU.typeId, menuitem.Id)
#         if concurrentCommand:
#             threading.Thread(target=concurrentCommand,
#                              args=concurrentArgs,
#                              kwargs=concurrentKWargs).start()
#         G2frame.ProcessEvent(menuevent)
#     except Exception as msg:
#         asyncResults['error'] = f'Error from InvokeMenuCommand: {msg}'


# def InvokeMenuCommand(G2frame,menuitem):
#     '''Executes a menu found in :func:`GetMenuCommand`

#     :param wx.Frame G2frame: reference to main GSASII frame
#     :param wx.MenuItem menuitem: reference to a menu item

#     This is intended to be called from wx.CallAfter(),
#     so results are placed in dict asyncResults which
#     must be defined before this is called:
#      * asyncResults['error'] has error information, if an error occurs

#     '''
#     asyncResults.clear()
#     try:
#         menuevent = wx.CommandEvent(wx.EVT_MENU.typeId, menuitem.Id)
#         G2frame.ProcessEvent(menuevent)
#     except Exception as msg:
#         asyncResults['error'] = f'Error from InvokeMenuCommand: {msg}'

# def getButtons(winname, btnList):
#     '''Get one or more button from a window name and a 
#     list of button names

#     This should be called by CallAfter
#     '''
#     def getAllChildren(parent):
#         '''Recursively finds all child windows of a given parent window.

#         :Returns: A list containing all descendant wx.Window objects.
#         '''
#         all_children = []
#         for child in parent.GetChildren():
#             all_children.append(child)
#             all_children.extend(getAllChildren(child))
#         return all_children
    
#     asyncResults.clear()
#     # find window
#     for w in wx.GetTopLevelWindows():
#         #print(w.GetTitle())
#         if w.GetTitle() == winname:
#             #print(f"Window {winname!r} found")
#             win = w
#             asyncResults['frame'] = win
#             break
#     else:
#         asyncResults['error'] = f"Window {winname!r} not found"        
#         return
#     # find each button
#     btns = []
#     for btnname in btnList:
#         for b in getAllChildren(win):
#             if not hasattr(b,'GetLabelText'): continue # not button
#             if not issubclass(b.__class__,wx.Button): continue # not button
#             if b.GetLabelText() == btnname:
#                 btns.append(b)
#                 break
#         else:
#             asyncResults['error'] = f'Button {btnname} not found'
#             return
#     asyncResults['buttons'] = btns

# def invokeButton(frame, btn):
#     '''"press" button from a window name and a button name

#     This should be called by CallAfter
#     '''
#     asyncResults.clear()
#     print(f'Invoking button {btn} id={btn.GetId()} lbl={btn.GetLabelText()}')
#     # why does this not work for "Set All" button?
#     #buttonevt = wx.PyCommandEvent(wx.EVT_BUTTON.typeId, btn.GetId())
#     #btn.GetTopLevelParent().ProcessEvent(buttonevt)
#     #frame.ProcessEvent(buttonevt)
#     #simulate mouse press on button
#     # G2frame = wx.GetApp().GetMainTopWindow()
#     # print('btn.GetPosition()',btn.GetPosition(),
#     #         '\nbtn.GetSize()',btn.GetSize(),
#     #         '\nframe.GetPosition()',frame.GetPosition(),
#     #         '\nG2frame.GetPosition()',G2frame.GetPosition(),)
#     # #pos = btn.GetPosition() + btn.GetSize()/2
#     # winpos = frame.GetPosition()
#     pos = frame.ClientToScreen(btn.GetPosition() + btn.GetSize()/2)
#     sim = wx.UIActionSimulator()
#     sim.MouseMove(pos.x,pos.y)
#     pause(0.05)
#     sim.MouseClick(wx.MOUSE_BTN_LEFT)
#     print('invoke done @',pos)
#     asyncResults['done'] = True

# def doInWindow(winname,bnlList):
#     pause(0.5) # wait for window to open
#     wx.CallAfter(getButtons,winname,bnlList)
#     pause(0.51)
#     print('getButtons',asyncResults)
#     if 'error' in asyncResults: return False
#     buttons,frame = asyncResults['buttons'],asyncResults['frame']
#     for btxt,b in zip(bnlList,buttons):
#         wx.CallAfter(invokeButton,frame,b)
#         count = 0
#         while len(asyncResults) == 0:
#             count += 1
#             pause(0.05)
#         print('invoke done',btxt,b,count)
#         #if 'error' in asyncResults: return False


        
#         #wx.CallAfter(InvokeMenuCommand,G2frame,menuitem,
#         #                 concurrentCommands,(G2frame,))
#         pause(0.1)
#         print(asyncResults)
#         if 'error' in asyncResults: return False
#     return
    
    # if False:
    #     # select tab on data window
    #     for i in range(G2frame.phaseDisplay.GetPageCount()):
    #         if G2frame.phaseDisplay.GetPageText(i) == 'Draw Atoms':
    #             G2frame.phaseDisplay.SetSelection(i)
    #             break
    #     else:
    #         print('Draw Atoms not found')
    #         return
    #     pause(0.1)

    #     # select menu command
    #     for i in range(G2frame.GetMenuBar().GetMenuCount()):
    #         if G2frame.GetMenuBar().GetMenuLabelText(i) == 'Edit Figure':
    #             menu = G2frame.GetMenuBar().GetMenu(i)
    #             break
    #     else:
    #         print('Edit Figure not found')
    #         return
    #     for i in menu.GetMenuItems():
    #         if menu.GetLabelText(i.Id) == 'Fill unit cell':
    #             menuitem = i
    #             break
    #     else:
    #         print('Fill unit cell not found')
    #         return
            
    #     #wx.CallLater(200,concurrentCommands,G2frame)
    #     #wx.CallAfter(concurrentCommands,G2frame)
    #     threading.Thread(target=concurrentCommands, args=(G2frame,)).start()
    #     G2frame.ProcessEvent(event)
    #     print('pause after menu command')
    #     pause(0.25)

    #     #breakpoint()

    #     #break
    #     #win.GetPosition()

    #     # go on to next phase entry
# def concurrentCommands(G2frame):
#     'finds a window and presses buttons in that window'
#     print('start concurrentCommands')
#     pause(0.5)
#     # button press find window named 'Select atoms for action'
#     # ['Select atoms for action', wx.Point(54, 191)]
#     winname = 'Select atoms for action'
#     pos = wx.Point(54, 191)
#     OKpos = wx.Point(175, 267)
#     for w in wx.GetTopLevelWindows():
#         if w.GetTitle() == winname:
#             print(f"Window {winname!r} found")
#             win = w
#             break
#     else:
#         print(f"Window {winname!r} not found")
#         return
#     winpos = win.GetPosition()
#     G2pos = G2frame.GetPosition()
#     print('modal,main',winpos,G2pos)
#     sim = wx.UIActionSimulator()
#     wx.CallAfter(sim.MouseMove,winpos.x+pos.x,winpos.y+pos.y)
#     pause(0.5)
#     print('before mouse click')
#     wx.CallAfter(sim.MouseClick,wx.MOUSE_BTN_LEFT)
#     pause(0.5)
#     print('mouse click')
#     wx.CallAfter(sim.MouseMove,winpos.x+OKpos.x,winpos.y+OKpos.y)
#     pause(0.5)
#     print('before mouse click')
#     wx.CallAfter(sim.MouseClick,wx.MOUSE_BTN_LEFT)
#     print('mouse click')

#     print('done concurrentCommands')
#     return

    
#     pause(0.5)
#     print('before mouse click')
#     wx.CallAfter(sim.MouseClick,wx.MOUSE_BTN_LEFT)
#     pause(0.5)
#     print('mouse click')


    # from importlib import reload
    # import playback
    # from playback import doAfterFileLoad
    # reload(playback)
    # from playback import doAfterFileLoad
    # import threading
    # #wx.CallAfter(doAfterFileLoad,G2frame)   # logging playback test
    # threading.Thread(target=doAfterFileLoad, args=(G2frame,)).start()
    # #doAfterFileLoad(G2frame)

#     def TestStuff(event): # playback
#         #============================================================
#         # animation code: open menu and invoke mouse command
#         from time import sleep
#         menupos = (132, 15)
#         cmdpos= (159, 186)
#         sim = wx.UIActionSimulator()
#         print('move to menu')
#         sim.MouseMove(menupos)
#         wx.Yield()
#         if sys.platform == "darwin":
#             dlg = wx.PopupTransientWindow(G2frame, wx.BORDER_RAISED)
#             wx.StaticText(dlg, label='Click on the "Save Project" entry')
#             dlg.SetPosition((cmdpos[0]+250,cmdpos[1]))
#             dlg.SetSize((250,40))
#             dlg.SetBackgroundColour(wx.YELLOW)
#             dlg.Show()
#             wx.Yield()
#             #breakpoint()
#             print('click on menu')
#             sim.MouseClick(wx.MOUSE_BTN_LEFT)
#             print('menu open')
#             wx.Yield()
#             dlg.Destroy()
#             return
#         # this is untested but might work on Linux/Windows
#         print('click on menu')
#         sim.MouseClick(wx.MOUSE_BTN_LEFT)
#         print('menu open')
#         for i in range(10): # move mouse slowly
#             nextpos = (np.array(cmdpos)-menupos)*(1+i)/10+menupos
#             print('move to ',nextpos)
#             sim.MouseMove(nextpos)
#             wx.Yield()
#             sleep(0.2)
#         # invoke the menu
#         sim.MouseClick(wx.MOUSE_BTN_LEFT)
#         return
#         #============================================================
#         # put text into a widget
#         typ = 'wxTextCtrl'
#         pos = (562, 69)
#         txt = "0.1234"
#         # get all windows in a frame
#         def get_all_windows(parent_window):
#             """
#             Recursively gets all child windows of a given wx.Window instance.

#             Args:
#             parent_window (wx.Window): The starting window to begin the traversal.

#             Returns:
#             list: A list containing all child windows found.
#             """
#             all_children = []
#             children = parent_window.GetChildren()
#             if children:
#                 for child in children:
#                     all_children.append(child)
#                     # Recursively call for each child
#                     all_children.extend(get_all_windows(child))
#             return all_children
#         #============================================================
#         # select a tree item
#         # need to find item
#         #item = GetGPXtreeItemId(G2frame,G2frame.root,'Phases')
#         #item = GetGPXtreeItemId(G2frame,item,'silicon')
#         item = GetGPXtreeItemId(G2frame,G2frame.root,'Notebook')
#         G2frame.GPXtree.SelectItem(item)
#         return
#         #============================================================
#         # create a lookup table for events once
#         # global G2EventLogger_evtCode
#         # if G2EventLogger_evtCode is None:
#         #     print('*** populating event table')
#         #     names = [n for n in dir(wx) if n.startswith('EVT_')]
#         #     evtdict = {n:getattr(wx,n) for n in names}
#         #     G2EventLogger_evtCode = {evt.typeId:n for n,evt in evtdict.items() if isinstance(evt, wx.PyEventBinder)}
#         # names = [n for n in dir(wx) if n.startswith('EVT_')]
#         # evtdict = {n:getattr(wx,n) for n in names}
#         # evtCode = {}
#         # for n,evt in evtdict.items():
#         #     if not isinstance(evt, wx.PyEventBinder):
#         #         #print('skip',n,evt)
#         #         continue
#         #     i = evt.typeId
#         #     if i in evtCode:
#         #         print('duplicate',i,n,evtCode[i])
#         #     else:
#         #         evtCode[i] = n
#         return
#         #============================================================
#         # executes 'Save project' from File Menu
#         # need to find menu
#         fileMenu = G2frame.dataMenuBars[0].GetMenu(0)
#         # need to find menu item
#         menuId = None
#         for o in fileMenu.GetMenuItems():
#             if o.GetItemLabelText() == 'Save project':
#                 menuId = o.Id
#                 break
#         else:
#             print('Not found')
#         itemId = eval('wx.EVT_MENU').typeId
#         sim = wx.CommandEvent(itemId,menuId)
#         wx.PostEvent(G2frame,sim)
#         #============================================================
#         #breakpoint()
#     btn.Bind(wx.EVT_BUTTON,TestStuff)
#     G2frame.dataWindow.SetDataSize()
# #    G2frame.SendSizeEvent()
    
    # def doPlayBack(self):
    #     '''Playback logging commands (eventually)
    #     '''
    #     return
    #     # use with /Users/toby/Scratch/logging/Aminoffite__R100055-9__Powder__DIF_File__10335.gpx
    #     cmdlist = [
    #         ('Tree_selection:',['Controls']),
    #         ('Tree_selection:',['Aminoffite', 'Phases']),
    #         ('Tree_selection:',['Controls']),
    #         ('Tree_selection:',['Aminoffite_shifted', 'Phases']),
    #         ]

    #     wx.GetApp().Yield()
    #     print('entering doPlayBack')
    #     for cmd in cmdlist:
    #         if 'Tree_selection' in cmd[0]:
    #             print('\n*** selecting tree item',reversed(cmd[1]))
    #             Id = GetGPXtreeItemId(self,self.root, cmd[1][-1])
    #             if Id == 0:
    #                 print('Tree item not found',cmd[1][-1])
    #                 return
    #             if len(cmd[1]) == 2:
    #                 Id = GetGPXtreeItemId(self,Id, cmd[1][0])
    #                 if Id == 0:
    #                     print('Tree item not found',cmd[1][0])
    #                     return
    #             print('selecting',Id)
    #             self.GPXtree.SelectItem(Id)
    #         print('before wait')
    #         for i in range(50):
    #             wx.GetApp().Yield()
    #             time.sleep(0.1)
    #         print('after wait')
    #     print('leaving doPlayBack')
        
