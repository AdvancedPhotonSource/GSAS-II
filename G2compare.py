#!/usr/bin/env python
# -*- coding: utf-8 -*-
#GSAS-II Data/Model Comparison
########### SVN repository information ###################
# $Date: $
# $Author: toby $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
'''
'''

#TODO: in OnDataTreeSelChanged want to plot patterns
# Make PWDR unique (use histlist)
# add graphics
# implement project

import sys
import os
import platform
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    try:
        import _pickle as cPickle
    except:
        print('Warning: failed to import the optimized Py3 pickle (_pickle)')
        import pickle as cPickle

import wx
import numpy as np
import matplotlib as mpl
try:
    import OpenGL as ogl
except ImportError:
    pass
import scipy as sp

import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 4154 $")
import GSASIIfiles as G2fil
import GSASIIplot as G2plt
import GSASIIctrlGUI as G2G
import GSASIIobj as G2obj

__version__ = '0.0.1'

def cPickleLoad(fp):
    if '2' in platform.python_version_tuple()[0]:
        return cPickle.load(fp)
    else:
       return cPickle.load(fp,encoding='latin-1')
            
def main(application):
    '''Start up the GSAS-II GUI'''                        
    knownVersions = ['2.7','3.6','3.7','3.8']
    if platform.python_version()[:3] not in knownVersions: 
        dlg = wx.MessageDialog(None, 
                'GSAS-II requires Python 2.7.x or 3.6+\n Yours is '+sys.version.split()[0],
                'Python version error',  wx.OK)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()
        sys.exit()
            
    application.main = MakeTopWindow(None)  # application.main is the main wx.Frame
    application.SetTopWindow(application.main)
    # save the current package versions
    application.main.PackageVersions = G2fil.get_python_versions([wx, mpl, np, sp, ogl])
    try:
        application.SetAppDisplayName('GSAS-II Compare')
    except:
        pass
    #application.GetTopWindow().SendSizeEvent()
    application.GetTopWindow().Show(True)
    return application.GetTopWindow()
    
class MakeTopWindow(wx.Frame):
    '''Define the main frame and its associated menu items
    '''
    def __init__(self, parent):
        size = wx.Size(700,450)
        wx.Frame.__init__(self, name='dComp', parent=parent,
            size=size,style=wx.DEFAULT_FRAME_STYLE, title='GSAS-II data/model comparison')
        # plot window
        self.plotFrame = wx.Frame(None,-1,'dComp Plots',size=size,
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.G2plotNB = G2plt.G2PlotNoteBook(self.plotFrame,G2frame=self)
        self.plotFrame.Show()
        # menus
        Frame = self.GetTopLevelParent() # same as self
        Menu = wx.MenuBar()
        File = wx.Menu(title='')
        Menu.Append(menu=File, title='&File')
        item = File.Append(wx.ID_ANY,'&Import project...\tCtrl+O','Open a GSAS-II project file (*.gpx)')            
        self.Bind(wx.EVT_MENU, self.onLoadGPX, id=item.GetId())
#        item = File.Append(wx.ID_ANY,'&Import selected...','Open a GSAS-II project file (*.gpx)')
#        self.Bind(wx.EVT_MENU, self.onLoadSel, id=item.GetId())

        self.Mode = wx.Menu(title='')
        Menu.Append(menu=self.Mode, title='&Mode')
        self.wxID_Mode = {}
        for m in "Histogram","Phase","Project":
            i = self.wxID_Mode[m] = wx.NewId()
            item = self.Mode.AppendRadioItem(i,m,'Display {}s'.format(m))
            self.Bind(wx.EVT_MENU, self.onRefresh, id=item.GetId())
        item = self.Mode.Append(wx.ID_ANY,'Set histogram filter','Set a filter for histograms to display')
        self.Bind(wx.EVT_MENU, self.onHistFilter, id=item.GetId())
        
        Frame.SetMenuBar(Menu)
        # status bar
        self.Status = self.CreateStatusBar()
        self.Status.SetFieldsCount(2)
        # split the frame and add the tree
        self.mainPanel = wx.SplitterWindow(self, wx.ID_ANY, style=wx.SP_LIVE_UPDATE|wx.SP_3D)
        self.mainPanel.SetMinimumPaneSize(100)
        self.treePanel = wx.Panel(self.mainPanel, wx.ID_ANY,
            style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
        
#        self.dataWindow = G2DataWindow(self.mainPanel)
        self.dataWindow = wx.ScrolledWindow(self.mainPanel)
        dataSizer = wx.BoxSizer(wx.VERTICAL)
        self.dataWindow.SetSizer(dataSizer)
        self.mainPanel.SplitVertically(self.treePanel, self.dataWindow, 200)
        self.Status.SetStatusWidths([200,-1])   # make these match?
        
        treeSizer = wx.BoxSizer(wx.VERTICAL)
        self.treePanel.SetSizer(treeSizer)
        self.GPXtree = G2G.G2TreeCtrl(id=wx.ID_ANY,
            parent=self.treePanel, size=self.treePanel.GetClientSize(),style=wx.TR_DEFAULT_STYLE )
        TreeId = self.GPXtree.Id

        treeSizer.Add(self.GPXtree,1,wx.EXPAND|wx.ALL,0)
        #self.GPXtree.Bind(wx.EVT_TREE_SEL_CHANGED,self.OnDataTreeSelChanged)
        self.GPXtree.Bind(wx.EVT_TREE_SEL_CHANGED,self.OnDataTreeSelChanged)
        # self.GPXtree.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK,self.OnDataTreeSelChanged)
        # self.GPXtree.Bind(wx.EVT_TREE_ITEM_COLLAPSED,
        #     self.OnGPXtreeItemCollapsed, id=TreeId)
        #self.GPXtree.Bind(wx.EVT_TREE_ITEM_EXPANDED,
        #     self.OnGPXtreeItemExpanded, id=TreeId)
        # self.GPXtree.Bind(wx.EVT_TREE_DELETE_ITEM,
        #     self.OnGPXtreeItemDelete, id=TreeId)
        # self.GPXtree.Bind(wx.EVT_TREE_KEY_DOWN,
        #     self.OnGPXtreeKeyDown, id=TreeId)
        # self.GPXtree.Bind(wx.EVT_TREE_BEGIN_RDRAG,
        #     self.OnGPXtreeBeginRDrag, id=TreeId)        
        # self.GPXtree.Bind(wx.EVT_TREE_END_DRAG,
        #     self.OnGPXtreeEndDrag, id=TreeId)        
        self.root = self.GPXtree.root        
        self.Bind(wx.EVT_CLOSE, lambda event: sys.exit())

        self.fileList = []  # list of files read for use in Reload
        self.histList = []  # list of histograms loaded for unique naming

        self.PWDRfilter = None
        
    def SelectGPX(self):
        '''Select a .GPX file to be read
        '''
        dlg = wx.FileDialog(self, 'Choose GSAS-II project file', 
                wildcard='GSAS-II project file (*.gpx)|*.gpx',style=wx.FD_OPEN)
        try:
            if dlg.ShowModal() != wx.ID_OK: return
            fil = os.path.splitext(dlg.GetPath())[0]+'.gpx'
        finally:
            dlg.Destroy()
        if os.path.exists(fil):
            self.fileList.append([fil,'GPX'])
            return fil
        else:
            print('File {} not found, skipping'.format(fil))
            return
        
    def getMode(self):
        '''returns the display mode (one of "Histogram","Phase","Project").
        Could return '?' in case of an error.
        '''
        for m in self.wxID_Mode:
            if self.Mode.FindItemById(self.wxID_Mode[m]).IsChecked():
                break
        else:
            m = '?'
        return m
    
    def onRefresh(self,event):
        '''reread all files, in response to a change in mode, etc.
        '''
        self.GPXtree.DeleteChildren(self.root)  # delete tree contents
        self.histList = []  # clear list of loaded histograms
        for fil,mode in self.fileList:
            self.loadFile(fil)
        
    def loadFile(self,fil):
        '''read or reread a file
        '''
        if self.getMode() == "Histogram":
            self.LoadPwdr(fil)
        elif self.getMode() == "Phase":
            self.LoadPhase(fil)
        elif self.getMode() == "Project":
            self.LoadProject(fil)
        else:
            print("mode not implemented")
            #raise Exception("mode not implemented")
        
    def onLoadGPX(self,event):
        '''Initial load of GPX file in response to a menu command
        '''
        fil = self.SelectGPX()
        if not fil: return
        if not os.path.exists(fil): return
        self.fileList.append([fil,'GPX'])
        self.loadFile(fil)

    def LoadPwdr(self,fil):
        '''Load PWDR entries from a .GPX file to the tree.
        see :func:`GSASIIIO.ProjFileOpen`
        '''
        G2frame = self
        filep = open(fil,'rb')
        shortname = os.path.splitext(os.path.split(fil)[1])[0]

        wx.BeginBusyCursor()
        histLoadList = []
        try:
            while True:
                try:
                    data = cPickleLoad(filep)
                except EOFError:
                    break
                if not data[0][0].startswith('PWDR'): continue
                if self.PWDRfilter is not None: # implement filter
                    if self.PWDRfilter not in data[0][0]: continue
                data[0][0] += ' ('
                data[0][0] += shortname
                data[0][0] += ')'
                histLoadList.append(data)
                        
        except Exception as errmsg:
            if GSASIIpath.GetConfigValue('debug'):
                print('\nError reading GPX file:',errmsg)
                import traceback
                print (traceback.format_exc())
            msg = wx.MessageDialog(G2frame,message="Error reading file "+
                str(fil)+". This is not a current GSAS-II .gpx file",
                caption="Load Error",style=wx.ICON_ERROR | wx.OK | wx.STAY_ON_TOP)
            msg.ShowModal()
        finally:
            filep.close()
            wx.EndBusyCursor()

        datum = None
        for i,data in enumerate(histLoadList):
            datum = data[0]
            datum[0] = G2obj.MakeUniqueLabel(datum[0],self.histList)
            Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=datum[0])
            self.histList.append(datum[0])
            # if 'ranId' not in datum[1][0]: # patch: add random Id if not present
            #     datum[1][0]['ranId'] = ran.randint(0,sys.maxsize)
            G2frame.GPXtree.SetItemPyData(Id,datum[1][:3])  #temp. trim off junk (patch?)
            for datus in data[1:]:
                sub = G2frame.GPXtree.AppendItem(Id,datus[0])
    #patch
                if datus[0] == 'Instrument Parameters' and len(datus[1]) == 1:
                    if datum[0].startswith('PWDR'):
                        datus[1] = [dict(zip(datus[1][3],zip(datus[1][0],datus[1][1],datus[1][2]))),{}]
                    else:
                        datus[1] = [dict(zip(datus[1][2],zip(datus[1][0],datus[1][1]))),{}]
                    for item in datus[1][0]:               #zip makes tuples - now make lists!
                        datus[1][0][item] = list(datus[1][0][item])
    #end patch
                G2frame.GPXtree.SetItemPyData(sub,datus[1])
        if datum: # was anything loaded?
            print('project load successful for {}'.format(datum[0]))
    #        G2frame.Status.SetStatusText('Mouse RB drag/drop to reorder',0)
    #    G2frame.SetTitleByGPX()
        self.GPXtree.Expand(self.root)
        
    def onHistFilter(self,event):
        'Load a filter string via a dialog in response to a menu event'
        lbl = ''
        if self.PWDRfilter is not None:
            lbl = self.PWDRfilter
        dlg = G2G.SingleStringDialog(self,'Set string',
                                'Set a string that must be in histogram name',
                                 lbl,size=(400,-1))
        if dlg.Show():
            if dlg.GetValue().strip() == '':
                self.PWDRfilter = None
            else:
                self.PWDRfilter = dlg.GetValue()
            dlg.Destroy()
            self.onRefresh(event)
        else:
            dlg.Destroy()

    def LoadPhase(self,fil):
        '''Load Phase entries from a .GPX file to the tree.
        see :func:`GSASIIIO.ProjFileOpen`
        '''
        G2frame = self
        filep = open(fil,'rb')
        shortname = os.path.splitext(os.path.split(fil)[1])[0]

        wx.BeginBusyCursor()
        Phases = None
        try:
            while True:
                try:
                    data = cPickleLoad(filep)
                except EOFError:
                    break
                if not data[0][0].startswith('Phase'): continue
                Phases = data
                #if self.PWDRfilter is not None: # implement filter
                #    if self.PWDRfilter not in data[0][0]: continue
                data[0][0] += ' ('
                if Phases:
                    data[0][0] += shortname
                    data[0][0] += ')'
                else: 
                    data[0][0] += shortname
                    data[0][0] += 'has no phases)'
                Phases = data
                break
                
        except Exception as errmsg:
            if GSASIIpath.GetConfigValue('debug'):
                print('\nError reading GPX file:',errmsg)
                import traceback
                print (traceback.format_exc())
            msg = wx.MessageDialog(G2frame,message="Error reading file "+
                str(fil)+". This is not a current GSAS-II .gpx file",
                caption="Load Error",style=wx.ICON_ERROR | wx.OK | wx.STAY_ON_TOP)
            msg.ShowModal()
        finally:
            filep.close()
            wx.EndBusyCursor()

        datum = None
        if Phases:
            datum = data[0]
            #datum[0] = G2obj.MakeUniqueLabel(datum[0],self.histList)
            Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=datum[0])
            G2frame.GPXtree.SetItemPyData(Id,datum[1])
            for datus in data[1:]:
                #datus[0] += ' ('
                #datus[0] += shortname
                #datus[0] += ')'
                sub = G2frame.GPXtree.AppendItem(Id,datus[0])
                G2frame.GPXtree.SetItemPyData(sub,datus[1])
        if datum: # was anything loaded?
            self.GPXtree.Expand(Id)
            print('project load successful for {}'.format(datum[0]))
    #        G2frame.Status.SetStatusText('Mouse RB drag/drop to reorder',0)
    #    G2frame.SetTitleByGPX()
        self.GPXtree.Expand(self.root)

    def LoadProject(self,fil):
        '''Load the Covariance entry from a .GPX file to the tree.
        see :func:`GSASIIIO.ProjFileOpen`
        '''
        G2frame = self
        filep = open(fil,'rb')
        shortname = os.path.splitext(os.path.split(fil)[1])[0]

        wx.BeginBusyCursor()
        Phases = None
        try:
            while True:
                try:
                    data = cPickleLoad(filep)
                except EOFError:
                    break
                if not data[0][0].startswith('Covariance'): continue
                Covar = data[0]
                GSASIIpath.IPyBreak_base()
                #if self.PWDRfilter is not None: # implement filter
                #    if self.PWDRfilter not in data[0][0]: continue
                Covar[0] = shortname + ' Covariance'
                Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=Covar[0])
                G2frame.GPXtree.SetItemPyData(Id,Covar[1])
                break
            else:
                print("{} does not have refinement results".format(shortname))
        except Exception as errmsg:
            if GSASIIpath.GetConfigValue('debug'):
                print('\nError reading GPX file:',errmsg)
                import traceback
                print (traceback.format_exc())
            msg = wx.MessageDialog(G2frame,message="Error reading file "+
                str(fil)+". This is not a current GSAS-II .gpx file",
                caption="Load Error",style=wx.ICON_ERROR | wx.OK | wx.STAY_ON_TOP)
            msg.ShowModal()
        finally:
            filep.close()
            wx.EndBusyCursor()
        self.GPXtree.Expand(self.root)
        
    def OnDataTreeSelChanged(self,event):
        item = event.GetItem()
        print('selected',item)
        
        #print(self.GPXtree._getTreeItemsList(item))
        # pltNum = self.G2plotNB.nb.GetSelection()
        # print(pltNum)
        # if pltNum >= 0:                         #to avoid the startup with no plot!
        #     self.G2plotNB.nb.GetPage(pltNum)
        #     NewPlot = False
        # else:
        #     NewPlot = True
        #if self.getMode() == "Histogram":
        #self.PatternId = self.PickId  = item
        #G2plt.PlotPatterns(self,plotType='PWDR',newPlot=NewPlot)
            
    # def OnGPXtreeItemExpanded(self,event):
    #     item = event.GetItem()
    #     print('expanded',item)
    #     print(self.GPXtree._getTreeItemsList(item))
    #     if item == self.root:
    #         event.StopPropagation()
    #     else:
    #         event.Skip(False)

if __name__ == '__main__':
    #if sys.platform == "darwin": 
    #    application = G2App(0) # create the GUI framework
    #else:
    application = wx.App(0) # create the GUI framework
    try:
        GSASIIpath.SetBinaryPath(True)
    except:
        print('Unable to run with current setup, do you want to update to the')
        try:
            if '2' in platform.python_version_tuple()[0]:            
                ans = raw_input("latest GSAS-II version? Update ([Yes]/no): ")
            else:
                ans = input("latest GSAS-II version? Update ([Yes]/no): ")                
        except:
            ans = 'no'
        if ans.strip().lower() == "no":
            import sys
            print('Exiting')
            sys.exit()
        print('Updating...')
        GSASIIpath.svnUpdateProcess()
    GSASIIpath.InvokeDebugOpts()
    Frame = main(application) # start the GUI
    argLoadlist = sys.argv[1:]
    if len(argLoadlist) == 0:
        argLoadlist = ['/Users/toby/Scratch/copy.gpx',
                       '/Users/toby/Scratch/CW_YAG.gpx']
    for arg in argLoadlist:
        fil = os.path.splitext(arg)[0] + '.gpx'
        if os.path.exists(fil):
            Frame.fileList.append([fil,'GPX'])
            Frame.loadFile(fil)
        else:
            print('File {} not found. Skipping'.format(fil))
    application.MainLoop()
