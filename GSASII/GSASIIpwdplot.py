# -*- coding: utf-8 -*-
'''
Classes and routines defined in :mod:`GSASIIpwdplot` follow. 
'''
from __future__ import division, print_function
import time
import copy
import math
import sys
import os.path
import numpy as np
import numpy.ma as ma
from . import GSASIIpath
# Don't depend on wx/matplotlib/scipy for scriptable; or for Sphinx docs
try:
    import wx
    import wx.aui
    import wx.glcanvas
except (ImportError, ValueError):
    print('GSASIIpwdplot: wx not imported')
try:
    import matplotlib as mpl
    if not mpl.get_backend():       #could be assigned by spyder debugger
        mpl.use('wxAgg')
    import matplotlib.figure as mplfig
except (ImportError, ValueError) as err:
    print('GSASIIpwdplot: matplotlib not imported')
    if GSASIIpath.GetConfigValue('debug'): print('error msg:',err)
from . import GSASIIdataGUI as G2gd
from . import GSASIIpwdGUI as G2pdG
from . import GSASIIlattice as G2lat
from . import GSASIImath as G2mth
from . import GSASIIctrlGUI as G2G
#import GSASIIobj as G2obj
from . import GSASIIplot as G2plt
import matplotlib.colors as mpcls
try:
    from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
except ImportError:
    from matplotlib.backends.backend_wx import FigureCanvas as Canvas
try:
    from matplotlib.backends.backend_agg import FigureCanvasAgg as hcCanvas
except ImportError:
    from matplotlib.backends.backend_agg import FigureCanvas as hcCanvas # standard name
except RuntimeError:  # happens during doc builds
    pass

# useful degree trig functions
sind = lambda x: math.sin(x*math.pi/180.)
cosd = lambda x: math.cos(x*math.pi/180.)
tand = lambda x: math.tan(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
acosd = lambda x: 180.*math.acos(x)/math.pi
atan2d = lambda x,y: 180.*math.atan2(y,x)/math.pi
atand = lambda x: 180.*math.atan(x)/math.pi
# numpy versions
npsind = lambda x: np.sin(x*np.pi/180.)
npcosd = lambda x: np.cos(x*np.pi/180.)
nptand = lambda x: np.tan(x*np.pi/180.)
npacosd = lambda x: 180.*np.arccos(x)/np.pi
npasind = lambda x: 180.*np.arcsin(x)/np.pi
npatand = lambda x: 180.*np.arctan(x)/np.pi
npatan2d = lambda x,y: 180.*np.arctan2(x,y)/np.pi
# misc global vars
Clip_on = GSASIIpath.GetConfigValue('Clip_on',True)
Gkchisq = chr(0x03C7)+chr(0xb2)
plotDebug = False
timeDebug = GSASIIpath.GetConfigValue('Show_timing',False)
obsInCaption = True # include the observed, calc,... items in the plot caption (PlotPatterns)
# options for publication-quality Rietveld plots
plotOpt = {}
plotOpt['labelSize'] = '11'
plotOpt['dpi'] = 600
plotOpt['width'] = 8.
plotOpt['height'] = 6.
plotOpt['Show'] = {}
plotOpt['legend'] = {}
plotOpt['colors'] = {}
plotOpt['format'] = None
plotOpt['initNeeded'] = True
plotOpt['lineList']  = ('obs','calc','bkg','zero','diff')
plotOpt['phaseLabels']  = {}
plotOpt['lineWid'] = '1'
plotOpt['saveCSV'] = False
plotOpt['CSVfile'] = None
for xy in 'x','y':
    for minmax in 'min','max':
        key = f'{xy}{minmax}'
        plotOpt[key] = 0.0
        plotOpt[key+'_use'] = False
partialOpts = {}  # options for display of partials
savedX = None     # contains the X values in plot units for "master" pattern, with mask

#### PlotPatterns ################################################################################
def ReplotPattern(G2frame,newPlot,plotType,PatternName=None,PickName=None):
    '''This does the same as PlotPatterns except that it expects the information
    to be plotted (pattern name, item picked in tree + eventually the reflection list)
    to be passed as names rather than references to wx tree items, defined as class entries
    '''
    if PatternName:
        pId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, PatternName)
        if pId:
            G2frame.PatternId = pId
        else:
            if GSASIIpath.GetConfigValue('debug'): print('PatternName not found',PatternName)
            return
    if PickName == PatternName:
        G2frame.PickId = G2frame.PatternId
    elif PickName:
        pId = G2gd.GetGPXtreeItemId(G2frame, G2frame.PatternId, PickName)
        if pId:
            G2frame.PickId = pId
        else:
            if GSASIIpath.GetConfigValue('debug'): print('PickName not found',PickName)
            return
    elif GSASIIpath.GetConfigValue('debug'):
        print('Possible PickId problem PickId=',G2frame.PickId)
    # for now I am not sure how to regenerate G2frame.HKL
    G2frame.HKL = []  # array of generated reflections
    G2frame.Extinct = [] # array of extinct reflections
    PlotPatterns(G2frame,plotType=plotType)

def plotVline(Page,Plot,Lines,Parms,pos,color,pick,style='dotted'):
    '''shortcut to plot vertical lines for limits & Laue satellites.
    Was used for extrapeaks'''
    if Page.plotStyle['qPlot']:
        Lines.append(Plot.axvline(2.*np.pi/G2lat.Pos2dsp(Parms,pos),color=color,
            picker=pick,pickradius=2.,linestyle=style))
    elif Page.plotStyle['dPlot']:
        Lines.append(Plot.axvline(G2lat.Pos2dsp(Parms,pos),color=color,
            picker=pick,pickradius=2.,linestyle=style))
    else:
        Lines.append(Plot.axvline(pos,color=color,
            picker=pick,pickradius=2.,linestyle=style))
        
def PlotPatterns(G2frame,newPlot=False,plotType='PWDR',data=None,
                     extraKeys=[],refineMode=False,indexFrom=''):
    '''Powder pattern plotting package - displays single or multiple powder 
    patterns as intensity vs 2-theta, q or TOF. Can display multiple patterns 
    as "waterfall plots" or contour plots. Log I plotting available.

    Note that information needed for plotting will be found in:

      * G2frame.PatternId: contains the tree item for the current histogram       
      * G2frame.PickId: contains the actual selected tree item (can be child 
        of histogram)
      * G2frame.HKL: used for tooltip display of hkl for a selected/generated 
        phase's reflections when mouse is moved to a reflection location; 
        HKL locations shown (usually as an orange line) in "Index Peak List"
        & "Unit Cells List" plots. 
        N.B. reflection tick markers are generated from each phase's 
        reflection list.
      * G2frame.Extinct: used for display of extinct reflections (in blue) 
        for generated reflections when "show extinct" is selected.
    '''
    global PlotList,IndxFrom
    IndxFrom = indexFrom
    def PublishPlot(event):
        msg = ""
        if 'PWDR' not in plottype:
            msg += " * only PWDR histograms can be used"
        if G2frame.Contour or not G2frame.SinglePlot:
            if msg: msg += '\n'
            msg += " * only when a single histogram is plotted"
        if Page.plotStyle['logPlot']:
            if msg: msg += '\n'
            msg += " * only when the intensity scale is linear/sqrt (not log)"
        if msg:
            msg = 'Publication export is only available under limited plot settings\n'+msg
            G2G.G2MessageBox(G2frame,msg,'Wrong plot settings')
            print(msg)
        elif G2frame.Weight:
            G2frame.Weight = False
            PlotPatterns(G2frame,plotType=plottype,extraKeys=extraKeys)
            PublishRietveldPlot(G2frame,Pattern,Plot,Page)
            G2frame.Weight = True
            PlotPatterns(G2frame,plotType=plottype,extraKeys=extraKeys)
            return
        else:
            PublishRietveldPlot(G2frame,Pattern,Plot,Page)

    def OnPlotKeyPress(event):
        try:        #one way to check if key stroke will work on plot
            Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
        except TypeError:
            G2frame.G2plotNB.status.SetStatusText('Select '+plottype+' pattern first',1)
            return
        newPlot = False
        if event.key == 'w':
            G2frame.Weight = not G2frame.Weight
            if not G2frame.Weight and not G2frame.Contour and 'PWDR' in plottype:
                G2frame.SinglePlot = True
            elif 'PWDR' in plottype: # Turning on Weight plot clears previous limits
                G2frame.FixedLimits['dylims'] = ['','']                
            newPlot = True
        elif event.key == 'X' and plottype == 'PWDR':
            G2frame.CumeChi = not G2frame.CumeChi 
        elif event.key == 'e' and plottype in ['SASD','REFD']:
            G2frame.ErrorBars = not G2frame.ErrorBars
        elif event.key == 'T' and 'PWDR' in plottype:
            Page.plotStyle['title'] = not Page.plotStyle.get('title',True)
        elif event.key == 'f' and 'PWDR' in plottype: # short or full length tick-marks
            Page.plotStyle['flTicks'] = not Page.plotStyle.get('flTicks',False)
        elif event.key == 'x'and 'PWDR' in plottype:
            Page.plotStyle['exclude'] = not Page.plotStyle['exclude']
        elif event.key == '.':
            Page.plotStyle['WgtDiagnostic'] = not Page.plotStyle.get('WgtDiagnostic',False)
            newPlot = True
        elif event.key == 'b' and plottype not in ['SASD','REFD'] and not Page.plotStyle['logPlot'] and not Page.plotStyle['sqrtPlot']:
            G2frame.SubBack = not G2frame.SubBack
        elif event.key == 'n':
            if G2frame.Contour:
                pass
            else:
                Page.plotStyle['logPlot'] = not Page.plotStyle['logPlot']
                if Page.plotStyle['logPlot']:
                    Page.plotStyle['sqrtPlot'] = False
                else:
                    Page.plotStyle['Offset'][0] = 0
                newPlot = True
        elif event.key == 's' and 'PWDR' in plottype:
            Page.plotStyle['sqrtPlot'] = not Page.plotStyle['sqrtPlot']
            if Page.plotStyle['sqrtPlot']:
                Page.plotStyle['logPlot'] = False
                G2frame.SubBack = False
            YmaxS = max(Pattern[1][1])
            if Page.plotStyle['sqrtPlot']:
                Page.plotStyle['delOffset'] = .02*np.sqrt(YmaxS)
                Page.plotStyle['refOffset'] = -0.1*np.sqrt(YmaxS)
                Page.plotStyle['refDelt'] = .1*np.sqrt(YmaxS)
            else:
                Page.plotStyle['delOffset'] = .02*YmaxS
                Page.plotStyle['refOffset'] = -0.1*YmaxS
                Page.plotStyle['refDelt'] = .1*YmaxS
            newPlot = True
        elif event.key == 'S' and 'PWDR' in plottype:
            choice = [m for m in mpl.cm.datad.keys()]+['GSPaired','GSPaired_r',]   # if not m.endswith("_r")
            choice.sort()
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Color scheme',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                G2frame.ContourColor = choice[sel]
            else:
                G2frame.ContourColor = GSASIIpath.GetConfigValue('Contour_color','Paired')
            dlg.Destroy()
            newPlot = True
        elif event.key == 'u' and (G2frame.Contour or not G2frame.SinglePlot):
            if G2frame.Contour:
                G2frame.Cmax = min(1.0,G2frame.Cmax*1.2)
            elif Page.plotStyle['Offset'][0] < 100.:
                Page.plotStyle['Offset'][0] += 1.
        elif event.key == 'd' and (G2frame.Contour or not G2frame.SinglePlot):
            if G2frame.Contour:
                G2frame.Cmax = max(0.0,G2frame.Cmax*0.8)
            elif Page.plotStyle['Offset'][0] > -100.:
                Page.plotStyle['Offset'][0] -= 1.
        elif (event.key == 'u' or event.key == 'd'
                  ) and G2frame.GPXtree.GetItemText(G2frame.PickId) in [
                'Index Peak List','Peak List']:
            # select the next or previous peak and highlight it
            if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Index Peak List':                grid = G2frame.indxPeaks
            elif G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Peak List':
                grid = G2frame.reflGrid
            selected = grid.GetSelectedRows()
            grid.ClearSelection()
            if event.key == 'd':
                if not selected: selected = [-1,]
                r = selected[0]+1
                if r >= grid.NumberRows: r = None
            else:
                if not selected: selected = [grid.NumberRows]
                r = selected[0]-1
                if r < 0: r = None
            if r is not None:
                grid.SelectRow(r,True)
                grid.MakeCellVisible(r,0)
            wx.CallAfter(grid.ForceRefresh)
        elif event.key == 'U':
            if G2frame.Contour:
                G2frame.Cmin += (G2frame.Cmax - G2frame.Cmin)/5.
            elif Page.plotStyle['Offset'][0] < 100.:
               Page.plotStyle['Offset'][0] += 10.
        elif event.key == 'D':
            if G2frame.Contour:
                G2frame.Cmin -= (G2frame.Cmax - G2frame.Cmin)/5.
            elif Page.plotStyle['Offset'][0] > -100.:
                Page.plotStyle['Offset'][0] -= 10.
        elif event.key == 'g':
            mpl.rcParams['axes.grid'] = not mpl.rcParams['axes.grid']
        elif event.key == 'l' and not G2frame.SinglePlot:
            Page.plotStyle['Offset'][1] -= 1.
        elif event.key == 'r' and not G2frame.SinglePlot:
            Page.plotStyle['Offset'][1] += 1.
        elif event.key == 'o':
            if G2frame.SinglePlot and not G2frame.Contour:
                global obsInCaption # include the observed, calc,... items in the plot caption (PlotPatterns)
                obsInCaption = not obsInCaption
            elif not G2frame.SinglePlot: 
                G2frame.Cmax = 1.0
                G2frame.Cmin = 0.0
                Page.plotStyle['Offset'] = [0,0]
        elif event.key == 'C' and 'PWDR' in plottype and G2frame.Contour:
            G2G.makeContourSliders(G2frame,Ymax,PlotPatterns,newPlot,plotType)
        elif event.key == 'c' and 'PWDR' in plottype:
            newPlot = True
            if not G2frame.Contour:
                G2frame.SinglePlot = False
                Page.plotStyle['Offset'] = [0.,0.]
                G2frame.FixedLimits['cylims'] = ['','']  # reset manual limits
            else:
                G2frame.SinglePlot = True                
            G2frame.Contour = not G2frame.Contour
        elif (event.key == 'p' and 'PWDR' in plottype and G2frame.SinglePlot):
            Page.plotStyle['partials'] = not Page.plotStyle['partials']
        elif (event.key == 'e' and 'PWDR' in plottype and G2frame.SinglePlot and ifLimits
                  and not G2frame.Contour):
            Page.excludeMode = not Page.excludeMode
            if Page.excludeMode:
                try: # fails from key menu
                    Page.startExclReg = event.xdata
                except AttributeError:
                    G2G.G2MessageBox(G2frame,'To create an excluded region, after clicking on "OK", move to the beginning of the region and press the "e" key. Then move to the end of the region and press "e" again','How to exclude')
                    return
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                axis.axvline(Page.startExclReg,color='b',dashes=(2,3))
                Page.canvas.draw() 
                Page.savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                y1, y2= Page.figure.axes[0].get_ylim()
                Page.vLine = Plot.axvline(Page.startExclReg,color='b',dashes=(2,3))
                Page.canvas.draw()
            else:
                Page.savedplot = None
                wx.CallAfter(PlotPatterns,G2frame,newPlot=False,
                    plotType=plottype,extraKeys=extraKeys)
                if abs(Page.startExclReg - event.xdata) < 0.1: return
                LimitId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits')
                data = G2frame.GPXtree.GetItemPyData(LimitId)
                mn = min(Page.startExclReg, event.xdata)
                mx = max(Page.startExclReg, event.xdata)
                data.append([mn,mx])
                G2pdG.UpdateLimitsGrid(G2frame,data,plottype)
            return
        elif event.key == 'a' and 'PWDR' in plottype and G2frame.SinglePlot and not (
                 Page.plotStyle['logPlot'] or Page.plotStyle['sqrtPlot'] or G2frame.Contour):
            # add a magnification region
            try:
                xpos = event.xdata
                if xpos is None: return  #avoid out of frame mouse position
                if 'Magnification' not in Pattern[0]:
                    Pattern[0]['Magnification'] = []
                try:
                    if Page.plotStyle['qPlot']:
                        xpos = G2lat.Dsp2pos(Parms,2.0*np.pi/xpos)
                    elif Page.plotStyle['dPlot']:
                        xpos = G2lat.Dsp2pos(Parms,xpos)
                except ValueError:
                    return
            except AttributeError: # invoked when this is called from dialog rather than key press
                xpos = (Pattern[1][0][-1]+Pattern[1][0][0])/2 # set to middle of pattern
            if not Pattern[0]['Magnification']:
                Pattern[0]['Magnification'] = [[None,1.]]
            Pattern[0]['Magnification'] += [[xpos,2.]]
            wx.CallAfter(G2gd.UpdatePWHKPlot,G2frame,plottype,G2frame.PatternId)
            return
        elif event.key == 'q' and not ifLimits: 
            newPlot = True
            if 'PWDR' in plottype:
                Page.plotStyle['qPlot'] = not Page.plotStyle['qPlot']
                Page.plotStyle['dPlot'] = False
                Page.plotStyle['chanPlot'] = False
            elif plottype in ['SASD','REFD']:
                Page.plotStyle['sqPlot'] = not Page.plotStyle['sqPlot']
        elif event.key == 'h' and G2frame.Contour:
            newPlot = True
            Page.plotStyle['qPlot'] = False
            Page.plotStyle['dPlot'] = False
            Page.plotStyle['chanPlot'] = not Page.plotStyle['chanPlot']
        elif event.key == 'e' and G2frame.Contour:
            newPlot = True
            G2frame.TforYaxis = not G2frame.TforYaxis
        elif event.key == 't' and 'PWDR' in plottype and not ifLimits:
            newPlot = True      
            Page.plotStyle['dPlot'] = not Page.plotStyle['dPlot']
            Page.plotStyle['qPlot'] = False
            Page.plotStyle['chanPlot'] = False
        elif event.key == 'm':
            if not G2frame.Contour:                
                G2frame.SinglePlot = not G2frame.SinglePlot                
            G2frame.Contour = False
            newPlot = True
        elif event.key == 'F' and not G2frame.SinglePlot:
            choices = G2gd.GetGPXtreeDataNames(G2frame,plotType)
            dlg = G2G.G2MultiChoiceDialog(G2frame,
                'Select dataset(s) to plot\n(select all or none to reset)', 
                'Multidata plot selection',choices)
            if dlg.ShowModal() == wx.ID_OK:
                G2frame.selections = []
                select = dlg.GetSelections()
                if select and len(select) != len(choices):
                    for Id in select:
                        G2frame.selections.append(choices[Id])
                else:
                    G2frame.selections = None
            dlg.Destroy()
            newPlot = True
        elif event.key in ['+','=']:
            G2frame.plusPlot = (G2frame.plusPlot+1)%3
        elif event.key == '/':
            Page.plotStyle['Normalize'] = not Page.plotStyle['Normalize']
            newPlot=True
        elif event.key == 'i' and G2frame.Contour:                  #for smoothing contour plot
            choice = ['nearest','bilinear','bicubic','spline16','spline36','hanning',
               'hamming','hermite','kaiser','quadric','catrom','gaussian','bessel',
               'mitchell','sinc','lanczos']
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Interpolation',choice)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                G2frame.Interpolate = choice[sel]
            else:
                G2frame.Interpolate = 'nearest'
            dlg.Destroy()
        elif event.key in [KeyItem[0] for KeyItem in extraKeys]:
            for KeyItem in extraKeys:
                if event.key == KeyItem[0]:
                    KeyItem[1]()
                    break
        elif event.key == 'v' and 'PWDR' in plottype and G2frame.SinglePlot:
            plotOpt['CSVfile'] = G2G.askSaveFile(G2frame,'','.csv',
                                        'Comma separated variable file')
            if plotOpt['CSVfile']: plotOpt['saveCSV'] = True
        else:
            #print('no binding for key',event.key)
            return
        wx.CallAfter(PlotPatterns,G2frame,newPlot=newPlot,plotType=plottype,extraKeys=extraKeys)
        
    def OnMotion(event):
        '''PlotPatterns: respond to motion of the cursor. The status line 
        will be updated with info based on the mouse position. 
        Also displays reflection labels as tooltips when mouse is over tickmarks
        '''
        global PlotList
        G2plt.SetCursor(Page)
        # excluded region animation
        if Page.excludeMode and Page.savedplot:
            if event.xdata is None or G2frame.GPXtree.GetItemText(
                    G2frame.GPXtree.GetSelection()) != 'Limits': # reset if out of bounds or not on limits
                Page.savedplot = None
                Page.excludeMode = False
                wx.CallAfter(PlotPatterns,G2frame,newPlot=False,
                    plotType=plottype,extraKeys=extraKeys)
                return
            else:
                Page.canvas.restore_region(Page.savedplot)
                Page.vLine.set_xdata([event.xdata,event.xdata])
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                axis.draw_artist(Page.vLine)
                Page.canvas.blit(axis.bbox)
                return
        elif Page.excludeMode or Page.savedplot: # reset if out of mode somehow
            Page.savedplot = None            
            Page.excludeMode = False
            wx.CallAfter(PlotPatterns,G2frame,newPlot=False,
                plotType=plottype,extraKeys=extraKeys)
            return
        if event.button and G2frame.Contour and G2frame.TforYaxis:
            ytics = imgAx.get_yticks()
            ytics = np.where(ytics<len(Temps),ytics,-1)
            imgAx.set_yticks(ytics)
            ylabs = [np.where(0<=i ,Temps[int(i)],' ') for i in ytics]
            imgAx.set_yticklabels(ylabs)            
        xpos = event.xdata
        if xpos is None: return  #avoid out of frame mouse position
        ypos = event.ydata
        try:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters')
            if not Id: return
            Parms,Parms2 = G2frame.GPXtree.GetItemPyData(Id)
            limx = Plot.get_xlim()
            dT = tolerance = np.fabs(limx[1]-limx[0])/100.
            if Page.plotStyle['qPlot'] and 'PWDR' in plottype:
                q = xpos
                if q <= 0:
                    G2frame.G2plotNB.status.SetStatusText('Q = %9.5f'%q)
                    return
                try:
                    dsp = 2.*np.pi/q
                    xpos = G2lat.Dsp2pos(Parms,2.0*np.pi/q)
                except ValueError:      #avoid bad value in asin beyond upper limit
                    G2frame.G2plotNB.status.SetStatusText('Q = %9.5f'%q)
                    return
                if 'T' in Parms['Type'][0]: # TOF
                    dT = Parms['difC'][1] * 2 * np.pi * tolerance / q**2
                elif 'E' in Parms['Type'][0]: # energy dispersive x-rays
                    pass    #for now
                else: # 'C' or  'B' in Parms['Type'][0] or 'PKS' in Parms['Type'][0]:
                    wave = G2mth.getWave(Parms)
                    dT = tolerance*wave*90./(np.pi**2*cosd(xpos/2))
            elif Page.plotStyle['chanPlot'] and G2frame.Contour:
                xpos = ma.getdata(X)[min(len(X)-1,int(xpos))]
                try:
                    dsp = G2lat.Pos2dsp(Parms,xpos)
                    q = 2.*np.pi/dsp
                except:
                    dsp = -1
                    q = -1
            elif plottype in ['SASD','REFD']:
                q = xpos
                if q <= 0:
                    G2frame.G2plotNB.status.SetStatusText('Q = %9.5f'%q)
                    return
                dsp = 2.*np.pi/q
            elif Page.plotStyle['dPlot']:
                dsp = xpos
                if dsp <= 0:
                    G2frame.G2plotNB.status.SetStatusText('d = %9.5f'%dsp)
                    return
                try:
                    q = 2.*np.pi/dsp
                    xpos = G2lat.Dsp2pos(Parms,dsp)
                except ValueError:      #avoid bad value
                    G2frame.G2plotNB.status.SetStatusText('d = %9.5f'%dsp)
                    return                
                dT = tolerance*xpos/dsp
            else:
                dsp = G2lat.Pos2dsp(Parms,xpos)
                q = 2.*np.pi/dsp
            statLine = ""
            if G2frame.Contour: #PWDR only
                try:
                    pNum = int(ypos+.5)
                    indx = abs(PlotList[pNum][1][0] - xpos).argmin() # closest point to xpos
                    val = 'int={:.3g}'.format(ma.getdata(PlotList[pNum][1][1])[indx])
                    if 'T' in Parms['Type'][0]:
                        statLine = 'TOF=%.3f d=%.5f Q=%.5f %s pattern ID=%d, %s'%(xpos,dsp,q,val,pNum,PlotList[pNum][-1])
                    else:
                        statLine = '2-theta=%.3f d=%.5f Q=%.5f %s pattern ID=%d, %s'%(xpos,dsp,q,val,pNum,PlotList[pNum][-1])
                except IndexError:
                    pass
            else:
                if 'T' in Parms['Type'][0]:
                    if Page.plotStyle['sqrtPlot']:
                        statLine = 'TOF = %9.3f d=%9.5f Q=%9.5f sqrt(Intensity) =%9.2f'%(xpos,dsp,q,ypos)
                    else:
                        statLine = 'TOF =%9.3f d=%9.5f Q=%9.5f Intensity =%9.2f'%(xpos,dsp,q,ypos)
                elif 'E' in Parms['Type'][0]:
                    statLine = 'Energy =%9.3f d=%9.5f Q=%9.5f sqrt(Intensity) =%9.2f'%(xpos,dsp,q,ypos)
                else:
                    if 'PWDR' in plottype:
                        ytmp = ypos
                        if Page.plotStyle['sqrtPlot']:
                            ytmp = ypos**2
                        statLine = '2-theta=%.3f d=%.5f Q=%.4f Intensity=%.2f'%(xpos,dsp,q,ytmp)
                    elif plottype == 'SASD':
                        statLine = 'q =%12.5g Intensity =%12.5g d =%9.1f'%(q,ypos,dsp)
                    elif plottype == 'REFD':
                        statLine = 'q =%12.5g Reflectivity =%12.5g d =%9.1f'%(q,ypos,dsp)
            zoomstat = Page.toolbar.get_zoompan()
            if zoomstat:
                statLine = "[" + zoomstat + "] " + statLine
            G2frame.G2plotNB.status.SetStatusText(statLine + IndxFrom,1)
            s = ''
            if G2frame.PickId:
                pickIdText = G2frame.GPXtree.GetItemText(G2frame.PickId)
            else:
                pickIdText = '?' # unexpected
            tickMarkList = ['Index Peak List','Unit Cells List','Reflection Lists']
            if pickIdText in tickMarkList and len(G2frame.HKL):
                found = []
                indx = -1
                if pickIdText in ['Index Peak List','Unit Cells List',]:
                    indx = -2
                # finds reflections within 1% of plot range in units of plot
                findx = np.where(np.fabs(np.array(G2frame.HKL).T[indx]-xpos) < dT/2.)
                found = G2frame.HKL[findx]
                if len(G2frame.Extinct):
                    G2frame.Extinct = np.array(G2frame.Extinct)
                    f2 = G2frame.Extinct[np.where(np.fabs(G2frame.Extinct.T[indx]-xpos) < dT/2.)] 
                    found = np.concatenate((found,f2))
                if found.shape[0]:
                    if len(found[0]) > 6:   #SS reflections
                        fmt = "{:.0f},{:.0f},{:.0f},{:.0f}"
                        n = 4
                    else:
                        fmt = "{:.0f},{:.0f},{:.0f}"
                        n = 3
                    for i,hkl in enumerate(found):
                        if i >= 3:
                            s += '\n...'
                            break
                        if s: s += '\n'
                        s += fmt.format(*hkl[:n])
            elif G2frame.itemPicked: # not sure when this will happen
                s = '%9.5f'%(xpos)
            Page.SetToolTipString(s)

        except TypeError:
            G2frame.G2plotNB.status.SetStatusText('Select '+plottype+' pattern first',1)
                
    def OnPress(event): #ugh - this removes a matplotlib error for mouse clicks in log plots
        np.seterr(invalid='ignore')
                                                   
    def onMoveDiffCurve(event):
        '''Respond to a menu command to move the difference curve. 
        '''
        if not DifLine[0]:
            print('No difference curve!')
            return
        G2frame.itemPicked = DifLine[0]
        G2frame.G2plotNB.Parent.Raise()
        OnPickPwd(None)

    def onMoveTopTick(event):
        '''Respond to a menu command to move the tick locations. 
        '''
        if len(Page.phaseList) == 0:
            print("there are tick marks (no phases)")
            return
        G2frame.itemPicked = Page.tickDict[Page.phaseList[0]]
        G2frame.G2plotNB.Parent.Raise()
        OnPickPwd(None)
                
    def onMoveTickSpace(event):
        '''Respond to a menu command to move the tick spacing. 
        '''
        if len(Page.phaseList) == 0:
            print("there are tick marks (no phases)")
            return
        G2frame.itemPicked = Page.tickDict[Page.phaseList[-1]]
        G2frame.G2plotNB.Parent.Raise()
        OnPickPwd(None)
        
    def onMovePeak(event):
        reflGrid = G2frame.reflGrid
        selectedPeaks = list(set([row for row,col in reflGrid.GetSelectedCells()] +
            reflGrid.GetSelectedRows()))
        if len(selectedPeaks) != 1:
            tbl = reflGrid.GetTable().data
            choices = [f"{i[0]:.2f}" for i in tbl]
            dlg = G2G.G2SingleChoiceDialog(G2frame,'Select peak to move',
                'select peak',choices)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    selectedPeaks = [dlg.GetSelection()]
            finally:
                dlg.Destroy()
        if len(selectedPeaks) != 1: return
        G2frame.itemPicked = G2frame.Lines[selectedPeaks[0]+2] # 1st 2 lines are limits
        G2frame.G2plotNB.Parent.Raise()
        OnPickPwd(None)

    def OnPickPwd(event):
        '''Respond to an item being picked. This usually means that the item
        will be dragged with the mouse, or sometimes the pick object is used 
        to create a peak or an excluded region
        '''
        def OnDragMarker(event):
            '''Respond to dragging of a plot Marker
            '''
            if event.xdata is None or event.ydata is None: return   # ignore if cursor out of window
            if G2frame.itemPicked is None: return # not sure why this happens, if it does
            Page.canvas.restore_region(savedplot)
            G2frame.itemPicked.set_data([event.xdata], [event.ydata])
            if G2frame.Weight:
                axis = Page.figure.axes[1]
            else:
                axis = Page.figure.gca()
            axis.draw_artist(G2frame.itemPicked)
            Page.canvas.blit(axis.bbox)
            
        def OnDragLine(event):
            '''Respond to dragging of a plot line
            '''
            if event.xdata is None: return   # ignore if cursor out of window
            if G2frame.itemPicked is None: return # not sure why this happens 
            Page.canvas.restore_region(savedplot)
            coords = G2frame.itemPicked.get_data()
            coords[0][0] = coords[0][1] = event.xdata
            coords = G2frame.itemPicked.set_data(coords)
            if G2frame.Weight:
                axis = Page.figure.axes[1]
            else:
                axis = Page.figure.gca()
            axis.draw_artist(G2frame.itemPicked)
            Page.canvas.blit(axis.bbox)
            
        def OnDragLabel(event):
            '''Respond to dragging of a HKL label
            '''
            if event.xdata is None: return   # ignore if cursor out of window
            if G2frame.itemPicked is None: return # not sure why this happens
            try:
                coords = list(G2frame.itemPicked.get_position())
                coords[1] = event.ydata
                data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
                data[0]['HKLmarkers'] = data[0].get('HKLmarkers',{})
                k1,k2 = G2frame.itemPicked.key
                if Page.plotStyle['sqrtPlot']:
                    data[0]['HKLmarkers'][k1][k2][0] = np.sign(coords[1])*coords[1]**2
                else:
                    data[0]['HKLmarkers'][k1][k2][0] = coords[1]
                Page.canvas.restore_region(savedplot)
                G2frame.itemPicked.set_position(coords)
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                axis.draw_artist(G2frame.itemPicked)
                Page.canvas.blit(axis.bbox)
            except:
                pass

        def OnDragTickmarks(event):
            '''Respond to dragging of the reflection tick marks
            '''
            if event.ydata is None: return   # ignore if cursor out of window
            if Page.tickDict is None: return # not sure why this happens, if it does
            Page.canvas.restore_region(savedplot)
            if Page.pickTicknum:
                refDelt = -(event.ydata-Page.plotStyle['refOffset'])/Page.pickTicknum
                refOffset = Page.plotStyle['refOffset']
            else:       #1st row of refl ticks
                refOffset = event.ydata
                refDelt = Page.plotStyle['refDelt']
            if G2frame.Weight:
                axis = Page.figure.axes[1]
            else:
                axis = Page.figure.gca()
            for pId,phase in enumerate(Page.phaseList):
                pos = refOffset - pId*refDelt
                coords = Page.tickDict[phase].get_data()
                coords[1][:] = pos
                Page.tickDict[phase].set_data(coords)
                axis.draw_artist(Page.tickDict[phase])
            Page.canvas.blit(axis.bbox)

        def OnDragDiffCurve(event):
            '''Respond to dragging of the difference curve. 
            '''
            if event.ydata is None: return   # ignore if cursor out of window
            if G2frame.itemPicked is None: return # not sure why this happens 
            Page.canvas.restore_region(savedplot)
            coords = G2frame.itemPicked.get_data()
            coords[1][:] += Page.diffOffset + event.ydata
            Page.diffOffset = -event.ydata
            G2frame.itemPicked.set_data(coords)
            Page.figure.gca().draw_artist(G2frame.itemPicked) #  Diff curve only found in 1-window plot
            Page.canvas.blit(Page.figure.gca().bbox)
            
        def DeleteHKLlabel(HKLmarkers,key):
            '''Delete an HKL label'''
            del HKLmarkers[key[0]][key[1]]
            G2frame.itemPicked = None
            PlotPatterns(G2frame,plotType=plottype,extraKeys=extraKeys)

        ####====== start of OnPickPwd
        plotNum = G2frame.G2plotNB.plotList.index('Powder Patterns')
        inXtraPeakMode = False   # Ignore if peak list menubar is not yet created
        try:
            inXtraPeakMode = G2frame.dataWindow.XtraPeakMode.IsChecked()
        except:
            pass
        global Page
        try:
            Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
        except TypeError:
            return
        if event is None: # called from a menu command rather than by click on mpl artist
            mouse = 1
            pick = G2frame.itemPicked
            ind = np.array([0])
        elif str(event.artist).startswith('Text'):  # respond to a right or left click on a HKL Label
            if G2frame.itemPicked is not None:  # only allow one selection 
                return
            pick = event.artist
            xpos,ypos = pick.get_position()
            if event.mouseevent.button == 3: # right click, delete HKL label
                # but only 1st of the picked items, if multiple
                G2frame.itemPicked = pick
                data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
                data[0]['HKLmarkers'] = data[0].get('HKLmarkers',{})
                # finish event processing before deleting the selection, so that
                # any other picked items are skipped
                wx.CallLater(100,DeleteHKLlabel,data[0]['HKLmarkers'],pick.key)
                return
            # prepare to drag HKL label (vertical position only)
            G2frame.itemPicked = pick
            pick.set_alpha(.3) # grey out text
            Page = G2frame.G2plotNB.nb.GetPage(plotNum)
            Page.canvas.draw() # refresh & save bitmap
            if G2frame.Weight:
                axis = Page.figure.axes[1]
            else:
                axis = Page.figure.gca()
            savedplot = Page.canvas.copy_from_bbox(axis.bbox)
            G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragLabel)
            pick.set_alpha(1.0)
            return
        else:
            if G2frame.itemPicked is not None: return
            pick = event.artist
            mouse = event.mouseevent
            xpos = pick.get_xdata()
            ypos = pick.get_ydata()
            ind = event.ind
            xy = list(list(zip(np.take(xpos,ind),np.take(ypos,ind)))[0])
            # convert from plot units
            xtick = xy[0] # selected tickmarck pos in 2theta/TOF or d-space (not Q)
            if Page.plotStyle['qPlot']:                              #qplot - convert back to 2-theta/TOF
                xy[0] = G2lat.Dsp2pos(Parms,2*np.pi/xy[0])
                xtick = xy[0]
            elif Page.plotStyle['dPlot']:                            #dplot - convert back to 2-theta/TOF
                xy[0] = G2lat.Dsp2pos(Parms,xy[0])
#            if Page.plotStyle['sqrtPlot']:
#                xy[1] = xy[1]**2
        if G2frame.PickId and G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Peak List':
            # Peak List: add peaks by clicking on points,
            # Or, by dragging: move peaks, limits, diff curve or (XtraPeaks only) tickmarks
            phname = ''
            try: 
                phname = str(pick).split('(',1)[1][:-1]
            except:
                pass
            if ind.all() != [0] and ObsLine[0].get_label() in str(pick):
                #picked a data point, add a new peak
                data = G2frame.GPXtree.GetItemPyData(G2frame.PickId)
                XY = G2mth.setPeakparms(Parms,Parms2,xy[0],xy[1],useFit=True)
                if inXtraPeakMode:
                    data['xtraPeaks'] = data.get('xtraPeaks',[])
                    data['xtraPeaks'].append(XY)
                    data['sigDict'] = {}    #now invalid
                else:
                    data['peaks'].append(XY)
                    data['sigDict'] = {}    #now invalid
                G2pdG.UpdatePeakGrid(G2frame,data)
                PlotPatterns(G2frame,plotType=plottype,extraKeys=extraKeys)
            elif inXtraPeakMode and phname in Page.phaseList:
                #picked a tick mark
                Page.pickTicknum = Page.phaseList.index(phname)
                resetlist = []
                for pId,phase in enumerate(Page.phaseList): # set the tickmarks to a lighter color
                    col = Page.tickDict[phase].get_color()
                    rgb = mpcls.ColorConverter().to_rgb(col)
                    rgb_light = [(2 + i)/3. for i in rgb]
                    resetlist.append((Page.tickDict[phase],rgb))
                    Page.tickDict[phase].set_color(rgb_light)
                    Page.tickDict[phase].set_zorder(99) # put on top
                Page.canvas.draw() # refresh with dimmed tickmarks 
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                for f,v in resetlist:  # reset colors back
                    f.set_zorder(0)
                    f.set_color(v) # reset colors back to original values
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragTickmarks)
                G2frame.itemPicked = pick
                return
            elif DifLine[0] is pick:
                # drag the difference curve
                G2frame.itemPicked = pick
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                Page.canvas.draw() # save bitmap
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                Page.diffOffset = Page.plotStyle['delOffset']
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragDiffCurve)
            else:
                #picked a peak list line, prepare to animate move of line
                G2frame.itemPicked = pick
                pick.set_linestyle(':') # set line as dotted
                Page = G2frame.G2plotNB.nb.GetPage(plotNum)
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                Page.canvas.draw() # refresh without dotted line & save bitmap
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragLine)
                pick.set_linestyle('--') # back to dashed
        elif G2frame.PickId and G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Limits':
            if len(event.artist.get_data()[0]) != 2 and not G2frame.ifSetLimitsMode: 
                return  # don't select a tickmark to drag. But in Limits/Excluded mode
                        # could select a point from the pattern here
            # Limits: add excluded region or move limits by use of menu command
            # and then pick a point
            # Or, drag line for limits/excluded region
            if ind.all() != [0]:                                    #picked a data point
                LimitId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits')
                limData = G2frame.GPXtree.GetItemPyData(LimitId)
                # Q & d not currently allowed on limits plot
                # if Page.plotStyle['qPlot']:                              #qplot - convert back to 2-theta
                #     xy[0] = G2lat.Dsp2pos(Parms,2*np.pi/xy[0])
                # elif Page.plotStyle['dPlot']:                            #dplot - convert back to 2-theta
                #     xy[0] = G2lat.Dsp2pos(Parms,xy[0])
                if G2frame.ifSetLimitsMode == 3:   # add an excluded region
                    excl = [0,0]
                    excl[0] = max(limData[1][0],min(xy[0],limData[1][1]))
                    excl[1] = excl[0]+0.1
                    limData.append(excl)
                elif G2frame.ifSetLimitsMode == 2: # set upper
                    limData[1][1] = max(xy[0],limData[1][0])
                elif G2frame.ifSetLimitsMode == 1:
                    limData[1][0] = min(xy[0],limData[1][1]) # set lower
                G2frame.ifSetLimitsMode = 0
                G2frame.CancelSetLimitsMode.Enable(False)
                G2frame.GPXtree.SetItemPyData(LimitId,limData)
                G2pdG.UpdateLimitsGrid(G2frame,limData,plottype)
                G2frame.GPXtree.SelectItem(LimitId)
                wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
                return
            else:                                         # picked a limit line
                # prepare to animate move of line
                G2frame.itemPicked = pick
                pick.set_linestyle(':') # set line as dotted
                Page = G2frame.G2plotNB.nb.GetPage(plotNum)
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                Page.canvas.draw() # refresh without dotted line & save bitmap
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragLine)
                pick.set_linestyle('--') # back to dashed
                
        elif G2frame.PickId and G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Unit Cells List':
            # By dragging lines: move limits
            if ind.all() == [0]:                         # picked a limit line
                # prepare to animate move of line
                G2frame.itemPicked = pick
                pick.set_linestyle(':') # set line as dotted
                Page = G2frame.G2plotNB.nb.GetPage(plotNum)
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                Page.canvas.draw() # refresh without dotted line & save bitmap
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragLine)
                pick.set_linestyle('--') # back to dashed

        elif G2frame.PickId and G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Models':
            if ind.all() != [0]:                                    #picked a data point
                LimitId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits')
                data = G2frame.GPXtree.GetItemPyData(LimitId)
                if mouse.button==1:
                    data[1][0] = min(xy[0],data[1][1])
                if mouse.button==3:
                    data[1][1] = max(xy[0],data[1][0])
                G2frame.GPXtree.SetItemPyData(LimitId,data)
                wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
            else:                                                   #picked a limit line
                G2frame.itemPicked = pick
        elif G2frame.PickId and (
                G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Reflection Lists'
                or 'PWDR' in G2frame.GPXtree.GetItemText(G2frame.PickId)
                ):
            G2frame.itemPicked = pick
            Page = G2frame.G2plotNB.nb.GetPage(plotNum)
            if G2frame.Weight:
                axis = Page.figure.axes[1]
            else:
                axis = Page.figure.gca()
            if DifLine[0] is G2frame.itemPicked:  # pick of difference curve
                Page.canvas.draw() # save bitmap
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                Page.diffOffset = Page.plotStyle['delOffset']
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragDiffCurve)
            elif G2frame.itemPicked in G2frame.MagLines: # drag of magnification marker
                pick.set_dashes((1,4)) # set line as dotted sparse
                Page.canvas.draw() # refresh without dotted line & save bitmap
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragLine)
                pick.set_dashes((1,1)) # back to dotted
            else:                         # pick of plot tick mark (is anything else possible?)
                #################################################################
                if event is not None and event.mouseevent.button == 3:
                    data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
                    font = int(data[0]['HKLconfig'].get('Font','8'))
                    angle = data[0]['HKLconfig'].get('Orientation',0)
                    G2frame.itemPicked = None # prevent tickmark repositioning
                    ph = event.artist.get_label()
                    if ph not in Phases:
                        print(f'Tickmark label ({ph}) not in Reflection Lists, very odd')
                        return
                    # find column for d-space or 2theta/TOF 
                    Super = 0
                    if Phases[ph].get('Super',False):
                        Super = 1
                    if Page.plotStyle['dPlot']:
                        indx = 4 + Super
                    else:
                        indx = 5 + Super
                    limx = Plot.get_xlim()
                    dT = tolerance = np.fabs(limx[1]-limx[0])/100.
                    if Page.plotStyle['qPlot']:
                        if 'T' in Parms['Type'][0]: # TOF
                            q = 2.*np.pi/G2lat.Pos2dsp(Parms,xpos)
                            dT = Parms['difC'][1] * 2 * np.pi * tolerance / q**2
                        elif 'E' in Parms['Type'][0]: # energy dispersive x-rays
                            pass    #for now
                        else: # 'C' or  'B' in Parms['Type'][0] or 'PKS' in Parms['Type'][0]:
                            wave = G2mth.getWave(Parms)
                            dT = tolerance*wave*90./(np.pi**2*cosd(xpos/2))
                    # locate picked tick markers
                    found_indices = np.where(np.fabs(Phases[ph]['RefList'][:,indx]-xtick) < dT/2.)
                    if len(found_indices) == 0: return
                    # get storage location for labeled reflections
                    data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
                    data[0]['HKLmarkers'] = data[0].get('HKLmarkers',{})
                    n = 3
                    if Super: n = 4
                    data[0]['HKLmarkers'][ph] = data[0]['HKLmarkers'].get(ph,{})
                    # position for 1st tick mark
                    ytick = ypos[0] - 3*(Plot.get_ylim()[1] - Plot.get_ylim()[0])/100
                    for i,refl in enumerate(Phases[ph]['RefList'][found_indices]):
                        key = f"{refl[5+Super]:.7g}"
                        data[0]['HKLmarkers'][ph][key] = data[0][
                            'HKLmarkers'][ph].get(key,[None,[]])
                        if Page.plotStyle['sqrtPlot']:
                            data[0]['HKLmarkers'][ph][key][0] = np.sign(ytick)*ytick**2 
                        else:
                            data[0]['HKLmarkers'][ph][key][0] = ytick
                        # add reflection if not already present
                        for hkl in data[0]['HKLmarkers'][ph][key][1]:
                            if sum(np.abs(refl[:n] - hkl)) == 0: break
                        else:
                            data[0]['HKLmarkers'][ph][key][1].append(tuple(np.rint(refl[:n]).astype(int)))
                        # offset next tick label position
                        ytick -= (1+angle)*3*(font/8.)*(Plot.get_ylim()[1] - Plot.get_ylim()[0])/100
                    wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
                    return
                #################################################################
                pick = str(G2frame.itemPicked).split('(',1)[1][:-1]
                if pick not in Page.phaseList: # picked something other than a tickmark
                    return
                Page.pickTicknum = Page.phaseList.index(pick)
                resetlist = []
                for pId,phase in enumerate(Page.phaseList): # set the tickmarks to a lighter color
                    col = Page.tickDict[phase].get_color()
                    rgb = mpcls.ColorConverter().to_rgb(col)
                    rgb_light = [(2 + i)/3. for i in rgb]
                    resetlist.append((Page.tickDict[phase],rgb))
                    Page.tickDict[phase].set_color(rgb_light)
                    Page.tickDict[phase].set_zorder(99) # put on top
                Page.canvas.draw() # refresh with dimmed tickmarks 
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                for f,v in resetlist:  # reset colors back
                    f.set_zorder(0)
                    f.set_color(v) # reset colors back to original values
                G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragTickmarks)
            
        elif G2frame.PickId and G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Background':
            # selected a fixed background point. Can move it or delete it.
            backPts = G2frame.dataWindow.wxID_BackPts
            for mode in backPts: # what menu is selected?
                if G2frame.dataWindow.BackMenu.FindItemById(backPts[mode]).IsChecked():
                    break
            # mode will be 'Add' or 'Move' or 'Del'
            if pick.get_marker() == 'D':
                # find the closest point
                backDict = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Background'))[1]
                d2 = [(x-xy[0])**2+(y-xy[1])**2 for x,y in backDict['FixedPoints']]
                G2frame.fixPtMarker = d2.index(min(d2))
                if mode == 'Move':
                    # animate move of FixedBkg marker
                    G2frame.itemPicked = pick
                    pick.set_marker('|') # change the point appearance
                    Page = G2frame.G2plotNB.nb.GetPage(plotNum)
                    if G2frame.Weight:
                        axis = Page.figure.axes[1]
                    else:
                        axis = Page.figure.gca()
                    Page.canvas.draw() # refresh with changed point & save bitmap
                    savedplot = Page.canvas.copy_from_bbox(axis.bbox)
                    G2frame.cid = Page.canvas.mpl_connect('motion_notify_event', OnDragMarker)
                    pick.set_marker('D') # put it back
                elif mode == 'Del':
                    del backDict['FixedPoints'][G2frame.fixPtMarker]
                    wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
                return
                
    def OnRelease(event):
        '''This is called when the mouse button is released when a plot object is dragged
        due to an item pick, or when invoked via a menu item (such as in onMoveDiffCurve),
        or for background points, which may be added/moved/deleted here.
        New peaks are also added here.
        '''
        plotNum = G2frame.G2plotNB.plotList.index('Powder Patterns')
        Page = G2frame.G2plotNB.nb.GetPage(plotNum)
        if G2frame.cid is not None:         # if there is a drag connection, delete it
            Page.canvas.mpl_disconnect(G2frame.cid)
            G2frame.cid = None
        if event.xdata is None or event.ydata is None: # ignore drag if cursor is outside of plot
#            if GSASIIpath.GetConfigValue('debug'): print('Ignoring drag, invalid pos:',event.xdata,event.ydata)
#            wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
            return
        if not G2frame.PickId:
            if GSASIIpath.GetConfigValue('debug'): print('Ignoring drag, G2frame.PickId is not set')
            return
        
        if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Background' and event.xdata:
            if Page.toolbar.AnyActive():    # prevent ops. if a toolbar zoom button pressed
                # after any mouse release event (could be a zoom), redraw magnification lines
                if magLineList: wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
                return 
            # Background page, deal with fixed background points
            if G2frame.SubBack or G2frame.Weight or G2frame.Contour or not G2frame.SinglePlot:
                return
            backDict = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Background'))[1]
            if 'FixedPoints' not in backDict: backDict['FixedPoints'] = []
            try:
                Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
            except TypeError:
                return
            # unit conversions
            xy = [event.xdata,event.ydata]
            try:
                if Page.plotStyle['qPlot']:                            #qplot - convert back to 2-theta
                    xy[0] = G2lat.Dsp2pos(Parms,2*np.pi/xy[0])
                elif Page.plotStyle['dPlot']:                          #dplot - convert back to 2-theta
                    xy[0] = G2lat.Dsp2pos(Parms,xy[0])
            except:
                return
            if Page.plotStyle['sqrtPlot']:
                xy[1] = xy[1]**2
            backPts = G2frame.dataWindow.wxID_BackPts
            for mode in backPts: # what menu is selected?
                if G2frame.dataWindow.BackMenu.FindItemById(backPts[mode]).IsChecked():
                    break
            if mode == 'Add':
                backDict['FixedPoints'].append(xy)
                if G2frame.Weight:
                    axis = Page.figure.axes[1]
                else:
                    axis = Page.figure.gca()
                axis.plot(event.xdata,event.ydata,'rD',clip_on=Clip_on,picker=True,pickradius=3.)
                Page.canvas.draw()
                return
            elif G2frame.itemPicked is not None: # end of drag in move
                backDict['FixedPoints'][G2frame.fixPtMarker] = xy
                G2frame.itemPicked = None
                wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
                return
        
        if G2frame.itemPicked is None:
            # after any mouse release event (could be a zoom) where nothing is selected,
            # redraw magnification lines
            if magLineList: wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
            return
        if DifLine[0] is G2frame.itemPicked:   # respond to dragging of the difference curve
            data = G2frame.GPXtree.GetItemPyData(G2frame.PickId)
            ypos = event.ydata
            Page.plotStyle['delOffset'] = -ypos
            G2frame.itemPicked = None
            wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
            return
        elif G2frame.itemPicked in G2frame.MagLines: # drag of magnification marker
            xpos = event.xdata
            try:
                if Page.plotStyle['qPlot']:                            #qplot - convert back to 2-theta
                    xpos = G2lat.Dsp2pos(Parms,2*np.pi/xpos)
                elif Page.plotStyle['dPlot']:                          #dplot - convert back to 2-theta
                    xpos = G2lat.Dsp2pos(Parms,xpos)
            except:
                return
            magIndex = G2frame.MagLines.index(G2frame.itemPicked)
            data = G2frame.GPXtree.GetItemPyData(G2frame.PickId)
            data[0]['Magnification'][magIndex][0] = xpos
            wx.CallAfter(G2gd.UpdatePWHKPlot,G2frame,plottype,G2frame.PatternId)
            return
        Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
        xpos = event.xdata
        phNum = None
        try:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Reflection Lists')
            pickPhase = str(G2frame.itemPicked).split('(',1)[1][:-1]
            phNum = Page.phaseList.index(pickPhase)
        except:
            pass
        if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Peak List'  and xpos and phNum is not None:
            # dragging tickmarks
            if phNum:
                Page.plotStyle['refDelt'] = -(event.ydata-Page.plotStyle['refOffset'])/phNum
            else:       #1st row of refl ticks
                Page.plotStyle['refOffset'] = event.ydata
        elif G2frame.GPXtree.GetItemText(G2frame.PickId) in ['Peak List','Limits','Unit Cells List'] and xpos:
            lines = []
            for line in G2frame.Lines: 
                lines.append(line.get_xdata()[0])
            try:
                lineNo = lines.index(G2frame.itemPicked.get_xdata()[0])
            except ValueError:
                lineNo = -1
            nxcl = len(exclLines)
            if  lineNo in [0,1] or lineNo in exclLines:
                LimitId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits')
                limits = G2frame.GPXtree.GetItemPyData(LimitId)
                Id = lineNo//2+1
                id2 = lineNo%2
                if Page.plotStyle['qPlot'] and 'PWDR' in plottype:
                    limits[Id][id2] = G2lat.Dsp2pos(Parms,2.*np.pi/xpos)
                elif Page.plotStyle['dPlot'] and 'PWDR' in plottype:
                    limits[Id][id2] = G2lat.Dsp2pos(Parms,xpos)
                else:
                    limits[Id][id2] = xpos
                if Id > 1 and limits[Id][0] > limits[Id][1]:
                        limits[Id].reverse()
                limits[1][0] = min(max(limits[0][0],limits[1][0]),limits[1][1])
                limits[1][1] = max(min(limits[0][1],limits[1][1]),limits[1][0])
                if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Limits':
                    G2pdG.UpdateLimitsGrid(G2frame,limits,plottype)
            elif lineNo > 1+nxcl:
                PeakId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List')
                peaks = G2frame.GPXtree.GetItemPyData(PeakId)
                inXtraPeakMode = False   # Ignore if peak list menubar is not yet created
                try:
                    inXtraPeakMode = G2frame.dataWindow.XtraPeakMode.IsChecked()
                except:
                    pass
                if inXtraPeakMode:
                    tbl = peaks['xtraPeaks']
                else:
                    tbl = peaks['peaks']
                if event.button == 3:
                    del tbl[lineNo-2-nxcl]
                else:
                    if Page.plotStyle['qPlot']:
                        tbl[lineNo-2-nxcl][0] = G2lat.Dsp2pos(Parms,2.*np.pi/xpos)
                    elif Page.plotStyle['dPlot']:
                        tbl[lineNo-2-nxcl][0] = G2lat.Dsp2pos(Parms,xpos)
                    else:
                        tbl[lineNo-2-nxcl][0] = xpos
                    peaks['sigDict'] = {}        #no longer valid
                G2pdG.UpdatePeakGrid(G2frame,peaks)
        elif G2frame.GPXtree.GetItemText(G2frame.PickId) in ['Models',] and xpos:
            lines = []
            for line in G2frame.Lines: 
                lines.append(line.get_xdata()[0])
            try:
                lineNo = lines.index(G2frame.itemPicked.get_xdata()[0])
            except ValueError:
                lineNo = -1
            if  lineNo in [0,1]:
                LimitId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits')
                data = G2frame.GPXtree.GetItemPyData(LimitId)
                data[1][lineNo] = xpos
                data[1][0] = min(max(data[0][0],data[1][0]),data[1][1])
                data[1][1] = max(min(data[0][1],data[1][1]),data[1][0])
        elif (G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Reflection Lists'
                or 'PWDR' in G2frame.GPXtree.GetItemText(G2frame.PickId)) and xpos:
            Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Reflection Lists')
            if Id:
                pick = str(G2frame.itemPicked).split('(',1)[1][:-1]
                if 'line' not in pick:       #avoid data points, etc.
                    if pick in Page.phaseList:
                        num = Page.phaseList.index(pick)
                        if num:
                            Page.plotStyle['refDelt'] = -(event.ydata-Page.plotStyle['refOffset'])/num
                        else:       #1st row of refl ticks
                            Page.plotStyle['refOffset'] = event.ydata
        elif GSASIIpath.GetConfigValue('debug'): # should not happen!
            print('How did we get here?')
            GSASIIpath.IPyBreak()
        wx.CallAfter(PlotPatterns,G2frame,plotType=plottype,extraKeys=extraKeys)
        G2frame.itemPicked = None

    def onSetPlotLim(event):
        '''Specify plot limits manually
        '''
        def onChecked(event):
            try:
                i = cbox.index(event.EventObject)
                showChecked(i)
            except:
                pass
        def showChecked(i):
            checked = cbox[i].GetValue()
            if not checked:
                # fake out validation to avoid ugly yellow 
                dbox[i].invalid = False
                dbox[i]._IndicateValidity()
            else: # reset validation
                dbox[i].ChangeValue(dbox[i].GetValue())
            dbox[i].Enable(checked)
        def applyLims(event):
            Page.toolbar.push_current()
            CurLims = {}
            CurLims['xlims'] = list(Plot.get_xlim())
            if G2frame.Weight:
                CurLims['ylims'] = list(Page.figure.axes[1].get_ylim())
                CurLims['dylims'] = list(Page.figure.axes[2].get_ylim())
            elif G2frame.Contour:
                CurLims['ylims'] = list(Plot.get_ylim())
                CurLims['cylims'] = list(Page.Img.get_clim())
            else:
                CurLims['ylims'] = list(Plot.get_ylim())
                CurLims['dylims'] = [0,0]
            for var in 'xlims','ylims','dylims','cylims':
                for i in range(2):
                    if not G2frame.UseLimits[var][i]: continue
                    try:
                        CurLims[var][i] = float(G2frame.FixedLimits[var][i])
                        CurLims[var][i] = float(G2frame.FixedLimits[var][i])
                    except:
                        pass
            Plot.set_xlim(CurLims['xlims'])
            if G2frame.Weight:
                Page.figure.axes[1].set_ylim(CurLims['ylims'])
                Page.figure.axes[2].set_ylim(CurLims['dylims'])
            elif G2frame.Contour:
                Plot.set_ylim(CurLims['ylims'])
                Page.Img.set_clim(CurLims['cylims'])
            else:
                Plot.set_ylim(CurLims['ylims'])
            Page.toolbar.push_current()
            Plot.figure.canvas.draw()
        
        # onSetPlotLim starts here
        dlg = wx.Dialog(G2frame.plotFrame,
                    style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(wx.StaticText(dlg,wx.ID_ANY,
                'Set Plot limits'
                ),0,wx.ALL)
        gsizer = wx.FlexGridSizer(cols=5,hgap=2,vgap=2)
        gsizer.Add(wx.StaticText(dlg,wx.ID_ANY,' '),0,wx.ALL)
        gsizer.Add(wx.StaticText(dlg,wx.ID_ANY,'use'),0,wx.ALL)
        gsizer.Add(wx.StaticText(dlg,wx.ID_ANY,'  min'),0,wx.ALL)
        gsizer.Add(wx.StaticText(dlg,wx.ID_ANY,'use'),0,wx.ALL)
        gsizer.Add(wx.StaticText(dlg,wx.ID_ANY,'  max'),0,wx.ALL)
        cbox = []
        dbox = []
        lblkeys = [(' x-axis ','xlims'),(' y-axis ','ylims')]
        if G2frame.Weight:
            lblkeys += [('(obs-calc)/sig ','dylims')]
        elif G2frame.Contour:
            lblkeys += [('contour','cylims')]

        for lbl,key in lblkeys:
            gsizer.Add(wx.StaticText(dlg,wx.ID_ANY,lbl),0,wx.ALL)
            for i in range(2):
                cbox.append(G2G.G2CheckBox(dlg,'',G2frame.UseLimits[key],i,
                                        OnChange=onChecked))
                dbox.append(G2G.ValidatedTxtCtrl(dlg,G2frame.FixedLimits[key],i,
                                                typeHint=float))
                gsizer.Add(cbox[-1])
                gsizer.Add(dbox[-1])
                showChecked(-1)
        vbox.Add(gsizer)
        vbox.Add((10,10),1,wx.ALL|wx.EXPAND,1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add((-1,-1),1,wx.ALL|wx.EXPAND,1)
        #btn = wx.Button(dlg, wx.ID_CLOSE)
        #btn.Bind(wx.EVT_BUTTON,lambda event:dlg.EndModal(wx.ID_CANCEL))
        #hbox.Add(btn)
        btn = wx.Button(dlg, wx.ID_ANY, label='Apply')
        btn.Bind(wx.EVT_BUTTON,applyLims)
        hbox.Add(btn)
        OKbtn = wx.Button(dlg, wx.ID_OK)
        OKbtn.Bind(wx.EVT_BUTTON,lambda event:dlg.EndModal(wx.ID_OK))
        hbox.Add(OKbtn)
        hbox.Add((-1,-1),1,wx.ALL|wx.EXPAND,1)
        vbox.Add(hbox,1,wx.ALL|wx.EXPAND,1)
        dlg.SetSizer(vbox)
        vbox.Fit(dlg)
        dlg.ShowModal()
        dlg.Destroy()
        applyLims(None) # apply limits
        
    def onPlotFormat(event):
        '''Change the appearance of the current plot'''
        changePlotSettings(G2frame,Plot)
        
    def refPlotUpdate(Histograms,cycle=None,restore=False):
        '''called to update an existing plot during a Rietveld fit; it only 
        updates the curves, not the reflection marks or the legend
        '''
        if restore:
            (G2frame.SinglePlot,G2frame.Contour,G2frame.Weight,
                G2frame.plusPlot,G2frame.SubBack,Page.plotStyle['logPlot']) = savedSettings
            return

        if plottingItem not in Histograms:
            histoList = [i for i in Histograms.keys() if i.startswith('PWDR ')]
            if len(histoList) == 0:
                print('Skipping plot, no PWDR item found!')
                return
            plotItem = histoList[0]
        else:
            plotItem = plottingItem
        xye = np.array(ma.getdata(Histograms[plotItem]['Data'])) # strips mask
        xye0 = Histograms[plotItem]['Data'][0]
        limits = Histograms[plotItem]['Limits']
        if Page.plotStyle['qPlot']:
            X = 2.*np.pi/G2lat.Pos2dsp(Parms,xye0)   # might want to consider caching this 
            Ibeg = np.searchsorted(X,2.*np.pi/G2lat.Pos2dsp(Parms,limits[1][0]))
            Ifin = np.searchsorted(X,2.*np.pi/G2lat.Pos2dsp(Parms,limits[1][1]))
        elif Page.plotStyle['dPlot']:
            X = G2lat.Pos2dsp(Parms,xye0)    # might want to consider caching this 
            Ibeg = np.searchsorted(X,G2lat.Pos2dsp(Parms,limits[1][1]))
            Ifin = np.searchsorted(X,G2lat.Pos2dsp(Parms,limits[1][0]))
        else:
            X = copy.deepcopy(xye0)
            Ibeg = np.searchsorted(X,limits[1][0])
            Ifin = np.searchsorted(X,limits[1][1])
        if Ibeg == Ifin: # if no points are within limits bad things happen 
            Ibeg,Ifin = 0,None
        if Page.plotStyle['sqrtPlot']:
            olderr = np.seterr(invalid='ignore') #get around sqrt(-ve) error
            Y = np.where(xye[1]>=0.,np.sqrt(xye[1]),-np.sqrt(-xye[1]))
            Z = np.where(xye[3]>=0.,np.sqrt(xye[3]),-np.sqrt(-xye[3]))
            W = np.where(xye[4]>=0.,np.sqrt(xye[4]),-np.sqrt(-xye[4]))
            #D = np.where(xye[5],(Y-Z),0.)-Page.plotStyle['delOffset']
            np.seterr(invalid=olderr['invalid'])
        else:
            Y = copy.copy(xye[1])
            Z = copy.copy(xye[3])
            W = copy.copy(xye[4])
            #D = xye[5]-Page.plotStyle['delOffset']  #powder background
        DZ = (xye[1]-xye[3])*np.sqrt(xye[2])
        DifLine[0].set_xdata(X[Ibeg:Ifin])
        DifLine[0].set_ydata(DZ[Ibeg:Ifin])
        try:
            lims = [min(DZ[Ibeg:Ifin]),max(DZ[Ibeg:Ifin])]
            if all(np.isfinite(lims)): Plot1.set_ylim(lims)
        except:
            pass
        CalcLine[0].set_xdata(X)
        ObsLine[0].set_xdata(X)
        BackLine[0].set_xdata(X)
        CalcLine[0].set_ydata(Z)
        ObsLine[0].set_ydata(Y)
        BackLine[0].set_ydata(W)
        if cycle:
            Title = '{} cycle #{}'.format(plotItem,cycle)
        else:
            Title = plotItem
        if Page.plotStyle.get('title',True):
            if Page.plotStyle['sqrtPlot']:
                Plot.set_title(r'$\sqrt{I}$ for '+Title)
            else:
                Plot.set_title(Title)
        Page.canvas.draw()
        
    def incCptn(string):
        '''Adds a underscore to "hide" a MPL object from the legend if 
        obsInCaption is False
        '''
        if obsInCaption:
            return string
        else:
            return '_'+string
    
    def Replot(*args):
        PlotPatterns(G2frame,plotType=plottype,extraKeys=extraKeys)
        
    def onHKLabelConfig(event):
        '''Configure all HKL markers in all phases'''
        def onDelAll(event):
            '''Delete all HKL markers in all phases'''
            data[0]['HKLmarkers'] = {}
            dlg.EndModal(wx.ID_OK)
            Replot()
        data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
        data[0]['HKLconfig'] = data[0].get('HKLconfig',{})
        dlg = wx.Dialog(G2frame,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(wx.StaticText(dlg,label='Configure Reflection (hkl) Labels'),
                      0,wx.ALIGN_CENTER_HORIZONTAL|wx.RIGHT|wx.LEFT|wx.BOTTOM,15)
        hsizer = wx.FlexGridSizer(cols=2,hgap=2,vgap=2)
        sizer.Add(hsizer)
        hsizer.Add(wx.StaticText(dlg,label='Label Orientation: '))
        data[0]['HKLconfig']['Orientation'] = data[0][
                'HKLconfig'].get('Orientation',0)
        choice = G2G.G2ChoiceButton(dlg,
                            ['Horizontal','Vertical'],
                            data[0]['HKLconfig'],'Orientation',
                            onChoice=Replot)
        hsizer.Add(choice)
        #
        hsizer.Add(wx.StaticText(dlg,label='Font size: '))
        data[0]['HKLconfig']['Font'] = data[0][
                'HKLconfig'].get('Font','8')
        choice = G2G.G2ChoiceButton(dlg,
                            ['6','8','10','12','14','16'],
                            None,None,
                            data[0]['HKLconfig'],'Font',
                            onChoice=Replot)
        hsizer.Add(choice)
        #
        hsizer.Add(wx.StaticText(dlg,label='Box transparency: '))
        data[0]['HKLconfig']['alpha'] = data[0][
                'HKLconfig'].get('alpha',2)
        choice = G2G.G2ChoiceButton(dlg,
                            ['0','25%','50%','75%','100%'],
                            data[0]['HKLconfig'],'alpha',
                            onChoice=Replot)
        hsizer.Add(choice)
        #
        btn = wx.Button(dlg, wx.ID_ANY,'Delete all labels')
        btn.Bind(wx.EVT_BUTTON, onDelAll)
        sizer.Add((-1,10))
        sizer.Add(btn, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5)
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(dlg, wx.ID_OK)
        btn.SetDefault()
        btn.Bind(wx.EVT_BUTTON, lambda x: dlg.EndModal(wx.ID_OK))
        btnsizer.AddButton(btn)
        btnsizer.Realize()
        sizer.Add(btnsizer, 0, wx.EXPAND|wx.ALL, 5)
        dlg.SetSizer(sizer)
        sizer.Fit(dlg)
        dlg.CenterOnParent()
        dlg.ShowModal()
        return
    
    def onPartialConfig(event):
        '''respond to Phase partials config menu command. 
        Computes partials if needed and then displays a window with options
        for each phase. Plot is refreshed each time a setting is changed.
        '''
        Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        phPartialFile = Controls['PhasePartials']
        if phPartialFile is None or not os.path.exists(phPartialFile):
            dlg = wx.MessageDialog(G2frame,
                'Phase partials (intensity profiles for each individual phase) have not been computed. Do you want to compute them?',
                    'Confirm compute',wx.YES_NO | wx.ICON_QUESTION)
            try:
                result = dlg.ShowModal()
            finally:
                dlg.Destroy()
            if result != wx.ID_YES: return
            G2frame.OnRefinePartials(None)
        if not Page.plotStyle['partials']:
            Page.plotStyle['partials'] = True
            Replot()
        configPartialDisplay(G2frame,Page.phaseColors,Replot)
            
    #### beginning PlotPatterns execution
    global exclLines,Page
    global DifLine # BHT: probably does not need to be global
    global Ymax
    global Pattern,mcolors,Plot,Page,imgAx,Temps
    global savedX
    plottype = plotType
    inXtraPeakMode = False   # Ignore if peak list menubar is not yet created
    try:
        inXtraPeakMode = G2frame.dataWindow.XtraPeakMode.IsChecked()
    except:
        pass

    # get powder pattern colors from config settings
    pwdrCol = {}
    for i in 'Obs_color','Calc_color','Diff_color','Bkg_color':
        pwdrCol[i] = '#' + GSASIIpath.GetConfigValue(i,getDefault=True)

    if not G2frame.PatternId:
        return
    if 'PKS' in plottype: # This is probably not used anymore; PlotPowderLines seems to be called directly
        G2plt.PlotPowderLines(G2frame)
        return
    if data is None:
        data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
    if G2frame.PickId and plottype not in ['SASD','REFD'] and 'PWDR' in G2frame.GPXtree.GetItemText(G2frame.PickId):
        publish = PublishPlot
    else:
        publish = None
    new,plotNum,Page,Plot,limits = G2frame.G2plotNB.FindPlotTab('Powder Patterns','mpl',publish=publish)
    if G2frame.ifSetLimitsMode and G2frame.GPXtree.GetItemText(G2frame.GPXtree.GetSelection()) == 'Limits':
        # note mode
        if G2frame.ifSetLimitsMode == 1:
            msg = 'Click on a point to define the location of the lower limit'
        elif G2frame.ifSetLimitsMode == 2:
            msg = 'Click on a point to define the location of the upper limit'
        elif G2frame.ifSetLimitsMode == 3:
            msg = 'Click on a point in the pattern to be excluded,\nthen drag or edit limits to adjust range'
        Page.figure.text(.02,.93, msg, fontsize=14, fontweight='bold')

    Page.excludeMode = False  # True when defining an excluded region
    Page.savedplot = None
#patch
    if 'Offset' not in Page.plotStyle and plotType in ['PWDR','SASD','REFD']:     #plot offset data
        Ymax = max(data[1][1])
        Page.plotStyle.update({'Offset':[0.0,0.0],'delOffset':0.02*Ymax,'refOffset':-0.1*Ymax,
            'refDelt':0.1*Ymax,})
#end patch
    if 'Normalize' not in Page.plotStyle:
        Page.plotStyle['Normalize'] = False
    # reset plot when changing between different data types
    try:
        G2frame.lastPlotType
    except:
        G2frame.lastPlotType = None
    if plotType == 'PWDR':
        try:
            Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
                G2frame.PatternId, 'Instrument Parameters'))
            if G2frame.lastPlotType != Parms['Type'][1]:
                if GSASIIpath.GetConfigValue('debug'): 
                    print('triggering newplot from G2frame.lastPlotType')
                Ymax = max(data[1][1])
                if Page.plotStyle['sqrtPlot']:
                    Page.plotStyle['delOffset'] = .02*np.sqrt(Ymax)
                    Page.plotStyle['refOffset'] = -0.1*np.sqrt(Ymax)
                    Page.plotStyle['refDelt'] = .1*np.sqrt(Ymax)
                else:
                    Page.plotStyle['delOffset'] = .02*Ymax
                    Page.plotStyle['refOffset'] = -0.1*Ymax
                    Page.plotStyle['refDelt'] = .1*Ymax
                newPlot = True
            G2frame.lastPlotType = Parms['Type'][1]
        except TypeError:       #bad id from GetGPXtreeItemId - skip
            pass
        
    try:
        G2frame.FixedLimits
    except:
        G2frame.FixedLimits = {i:['',''] for i in ('xlims','ylims','dylims','cylims')}
    try:
        G2frame.UseLimits
    except:
        G2frame.UseLimits = {i:[False,False] for i in ('xlims','ylims','dylims','cylims')}
    #=====================================================================================
    # code to setup for plotting Rietveld results. Turns off multiplot,
    # sqrtplot, turn on + and weight plot, but sqrtPlot qPlot and dPlot are not changed.
    # Magnification regions are ignored.
    # the last-plotted histogram (from G2frame.PatternId) is used for this plotting
    #    (except in seq. fitting)
    # Returns a pointer to refPlotUpdate, which is used to update the plot when this
    # returns
    if refineMode:
        plottingItem = G2frame.GPXtree.GetItemText(G2frame.PatternId)
        # save settings to be restored after refinement with repPlotUpdate({},restore=True)
        savedSettings = (G2frame.SinglePlot,G2frame.Contour,G2frame.Weight,
                            G2frame.plusPlot,G2frame.SubBack,Page.plotStyle['logPlot'])
        G2frame.SinglePlot = True
        G2frame.Contour = False
        G2frame.Weight = True
        G2frame.plusPlot = 1
        G2frame.SubBack = False
        Page.plotStyle['logPlot'] = False
        # is the selected histogram in the refinement? if not pick the 1st to show
        # instead select the first powder pattern and plot it
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        if plottingItem not in Histograms:  
            # current plotted item is not in refinement
            histoList = [i for i in Histograms.keys() if i.startswith('PWDR ')]
            if len(histoList) != 0:
                plottingItem = histoList[0]
                G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, plottingItem)
                data = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
                G2frame.GPXtree.SelectItem(G2frame.PatternId)
                PlotPatterns(G2frame,True,plotType,None,extraKeys)
    #=====================================================================================
    if not new:
        G2frame.xylim = copy.copy(limits)
    else:
        if plottype in ['SASD','REFD']:
            Page.plotStyle['logPlot'] = True
            G2frame.ErrorBars = True
        newPlot = True
        G2frame.Cmin = 0.0
        G2frame.Cmax = 1.0
        Page.canvas.mpl_connect('motion_notify_event', OnMotion)
        Page.canvas.mpl_connect('pick_event', OnPickPwd)
        Page.canvas.mpl_connect('button_release_event', OnRelease)
        Page.canvas.mpl_connect('button_press_event',OnPress)
        Page.bindings = []
    # redo OnPlotKeyPress binding each time the Plot is updated
    # since needs values that may have been changed after 1st call
    for b in Page.bindings:
        Page.canvas.mpl_disconnect(b)
    Page.bindings = []
    Page.bindings.append(Page.canvas.mpl_connect('key_press_event', OnPlotKeyPress))
    if not G2frame.PickId:
        print('No plot, G2frame.PickId,G2frame.PatternId=',G2frame.PickId,G2frame.PatternId)
        return
    elif 'PWDR' in G2frame.GPXtree.GetItemText(G2frame.PickId):
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        Phases = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Reflection Lists'))
        Page.phaseList = sorted(Phases.keys()) # define an order for phases (once!)
        G2frame.Bind(wx.EVT_MENU, onMoveDiffCurve, id=G2frame.dataWindow.moveDiffCurve.GetId())
        G2frame.Bind(wx.EVT_MENU, onMoveTopTick, id=G2frame.dataWindow.moveTickLoc.GetId())
        G2frame.Bind(wx.EVT_MENU, onMoveTickSpace, id=G2frame.dataWindow.moveTickSpc.GetId())
        G2frame.Bind(wx.EVT_MENU, onSetPlotLim, id=G2frame.dataWindow.setPlotLim.GetId())
        G2frame.Bind(wx.EVT_MENU, onPlotFormat, id=G2frame.dataWindow.setPlotFmt.GetId())
        G2frame.Bind(wx.EVT_MENU, onHKLabelConfig, id=G2G.wxID_CHHKLLBLS)
        if len(G2frame.Refine): # extend state for new menus to match main
            state = G2frame.Refine[0].IsEnabled()
        else:
            state = False
        G2frame.Bind(wx.EVT_MENU, onPartialConfig, id=G2G.wxID_CHPHPARTIAL)
        G2frame.PartialConfig.Enable(state)
        G2frame.Bind(wx.EVT_MENU, G2frame.OnSavePartials, id=G2G.wxID_PHPARTIALCSV)
        G2frame.PartialCSV.Enable(state) # disabled during sequential fits
        G2frame.dataWindow.moveDiffCurve.Enable(False)
        G2frame.dataWindow.moveTickLoc.Enable(False)
        G2frame.dataWindow.moveTickSpc.Enable(False)
    elif 'SAS' in plottype:
        Page.phaseList = Phases = []
    elif G2frame.GPXtree.GetItemText(G2frame.PickId) in ['Reflection Lists','Limits']:
        # get reflection positions for other places where they are displayed
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        Phases = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Reflection Lists'))
        Page.phaseList = sorted(Phases.keys()) # define an order for phases (once!)
    elif G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Peak List':
        G2frame.Bind(wx.EVT_MENU, onMovePeak, id=G2frame.dataWindow.movePeak.GetId())
        Page.phaseList = Phases = []
        if inXtraPeakMode:
            Phases = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId,'Reflection Lists'))
            Page.phaseList = sorted(Phases.keys()) # define an order for phases (once!)
    else:
        Page.phaseList = Phases = []
    # assemble a list of validated colors for tickmarks
    valid_colors = []
    invalid_colors = []
    for color in GSASIIpath.GetConfigValue('Ref_Colors',getDefault=True).split():
        try:
            if mpl.colors.is_color_like(color):
                valid_colors.append(color)
            else:
                invalid_colors.append(color)
        except:
            pass
    if invalid_colors and not hasattr(Page,'phaseColors'): # show error once
        print(f'**** bad color code(s): "{", ".join(invalid_colors)}" - redo Preferences/Ref_Colors ****')
    if len(valid_colors) < 3:
        refColors=['b','r','c','g','m','k']
    else:
        refColors = valid_colors
    if not hasattr(Page,'phaseColors'): Page.phaseColors = {}
    for i,p in enumerate(Phases):
        Page.phaseColors[p] = Page.phaseColors.get(p,refColors[i%len(refColors)])
    # save information needed to reload from tree and redraw
    try:
        if not refineMode:
            kwargs={'PatternName':G2frame.GPXtree.GetItemText(G2frame.PatternId)}
            if G2frame.PickId:
                kwargs['PickName'] = G2frame.GPXtree.GetItemText(G2frame.PickId)
            wx.CallAfter(G2frame.G2plotNB.RegisterRedrawRoutine(G2frame.G2plotNB.lastRaisedPlotTab,ReplotPattern,
                (G2frame,newPlot,plotType),kwargs))
    except:         #skip a C++ error
        pass
    # now start plotting
    G2frame.G2plotNB.status.DestroyChildren() #get rid of special stuff on status bar
    G2frame.G2plotNB.status.SetStatusText(IndxFrom,1)
    # TODO: figure out why the SetHelpButton creates a second tab line (BHT, Mac, wx4.1)
    #G2frame.G2plotNB.SetHelpButton(G2frame.dataWindow.helpKey)
    Page.tickDict = {}
    DifLine = ['']
    ifLimits = False
    if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Limits':
        ifLimits = True
        Page.plotStyle['qPlot'] = False
        Page.plotStyle['dPlot'] = False
    # keys in use for graphics control:
    #    a,b,c,d,e,f,g,i,l,m,n,o,p,q,r,s,t,u,w,x, (unused: j, k, y, z)
    #    also: +,/, C,D,S,U
    if G2frame.Contour:
        Page.Choice = (' key press','b: toggle subtract background',
            'd: lower contour max','u: raise contour max',
            'D: lower contour min','U: raise contour min',
            'o: reset contour limits','g: toggle grid',
            'i: interpolation method','S: color scheme','c: contour off',
            'e: toggle temperature for y-axis','s: toggle sqrt plot',
            'w: toggle w(Yo-Yc) contour plot','h: toggle channel # plot',
            'q: toggle Q plot','t: toggle d-spacing plot',
            'C: contour plot control window',
            )
    else:
        if 'PWDR' in plottype:
            Page.Choice = [' key press',
                'a: add magnification region','b: toggle subtract background',
                'c: contour on','x: toggle excluded regions','T: toggle plot title',
                'f: toggle full-length ticks','g: toggle grid',
                'X: toggle cumulative chi^2',
                'm: toggle multidata plot','n: toggle log(I)',]
            if obsInCaption:
                Page.Choice += ['o: remove obs, calc,... from legend',]
            else:
                Page.Choice += ['o: add obs, calc,... to legend',]
            if ifLimits:
                Page.Choice += ['e: create excluded region',
                        's: toggle sqrt plot','w: toggle (Io-Ic)/sig plot',
                        '+: toggle obs line plot']
            else:
                Page.Choice += [
                        'q: toggle Q plot','t: toggle d-spacing plot',
                        's: toggle sqrt plot','w: toggle (Io-Ic)/sig plot',
                        '+: toggle obs line plot']
            if G2frame.GPXtree.GetItemText(G2frame.PickId) in ['Peak List','Index Peak List']:
                Page.Choice += ['d: highlight next peak in list']
                Page.Choice += ['u: highlight previous peak in list']
            if Page.plotStyle['sqrtPlot'] or Page.plotStyle['logPlot']:
                del Page.Choice[1]
                del Page.Choice[1]
            elif not G2frame.SinglePlot:
                del Page.Choice[1]

            if not G2frame.SinglePlot:
                Page.Choice = Page.Choice+ \
                    ['u/U: offset up/10x','d/D: offset down/10x','l: offset left','r: offset right',
                     'o: reset offset','F: select data','/: normalize']
            else:
                Page.Choice = Page.Choice+ ['p: toggle partials (if available)',]
            if G2frame.SinglePlot:
                Page.Choice += ['v: CSV output of plot']
        elif plottype in ['SASD','REFD']:
            Page.Choice = [' key press',
                'b: toggle subtract background file','g: toggle grid',
                'm: toggle multidata plot','n: toggle semilog/loglog',
                'q: toggle S(q) plot','w: toggle (Io-Ic)/sig plot','+: toggle obs line plot',]
            if not G2frame.SinglePlot:
                Page.Choice = Page.Choice+ \
                    ['u: offset up','d: offset down','l: offset left',
                     'r: offset right','o: reset offset',]
                    
    for KeyItem in extraKeys:
        Page.Choice = Page.Choice + [KeyItem[0] + ': '+KeyItem[2],]
    magLineList = [] # null value indicates no magnification
    Page.toolbar.updateActions = None # no update actions
    G2frame.cid = None
    Page.keyPress = OnPlotKeyPress
    # assemble a list of validated colors (not currently needed)
    # valid_colors = []
    # invalid_colors = []
    # try:
    #     colors = GSASIIpath.GetConfigValue('Plot_Colors').split()
    #     for color in colors:
    #         #if color not in ['k','r','g','b','m','c']:
    #         if mpl.colors.is_color_like(color):
    #             valid_colors.append(color)
    #         else:
    #             invalid_colors.append(color)
    # except:
    #     pass
    # if invalid_colors:
    #     print(f'**** bad color code(s): "{", ".join(invalid_colors)}" - redo Preferences/Plot Colors ****')
    # if len(valid_colors) < 3:
    #     colors = ['b','g','r','c','m','k']
    # else:
    #     colors = valid_colors
    Lines = []
    exclLines = []
    time0 = time.time()
    if G2frame.SinglePlot and G2frame.PatternId:
        try:
            Pattern = G2frame.GPXtree.GetItemPyData(G2frame.PatternId)
            Pattern.append(G2frame.GPXtree.GetItemText(G2frame.PatternId))
            PlotList = [Pattern,]
            # PId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Background')
            # Pattern[0]['BackFile'] = ['',-1.0,False]
            # if PId:
            #     Pattern[0]['BackFile'] =  G2frame.GPXtree.GetItemPyData(PId)[1].get('background PWDR',['',-1.0,False])
            Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
                G2frame.PatternId, 'Instrument Parameters'))
            Sample = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Sample Parameters'))
            Limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))
            ParmList = [Parms,]
            SampleList = [Sample,]
            LimitsList = [Limits,]
            Title = data[0].get('histTitle')
            if not Title: 
                Title = Pattern[-1]
        except AttributeError:
            pass
    else:     #G2frame.selection   
        Title = os.path.split(G2frame.GSASprojectfile)[1]
        if G2frame.selections is None:
            choices = G2gd.GetGPXtreeDataNames(G2frame,plotType)
        else:
            choices = G2frame.selections
        PlotList = []
        ParmList = []
        SampleList = []
        LimitsList = []
        Temps = []
        # loop through tree looking for matching histograms to plot
        Id, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
        while Id:
            name = G2frame.GPXtree.GetItemText(Id)
            pid = Id
            Id, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
            if name not in choices: continue
            Pattern = G2frame.GPXtree.GetItemPyData(pid)
            if len(Pattern) < 3:                    # put name on end if needed
                Pattern.append(G2frame.GPXtree.GetItemText(pid))
            if 'Offset' not in Page.plotStyle:     #plot offset data
                Ymax = max(Pattern[1][1])
                Page.plotStyle.update({'Offset':[0.0,0.0],'delOffset':0.02*Ymax,'refOffset':-0.1*Ymax,'refDelt':0.1*Ymax,})
            # PId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Background')
            # Pattern[0]['BackFile'] = ['',-1.0,False]
            # if PId:
            #     Pattern[0]['BackFile'] =  G2frame.GPXtree.GetItemPyData(PId)[1].get('background PWDR',['',-1.0,False])
            PlotList.append(Pattern)
            ParmList.append(G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
                pid,'Instrument Parameters'))[0])
            SampleList.append(G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
                pid, 'Sample Parameters')))
            LimitsList.append(G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,
                pid, 'Limits')))
            Temps.append('%.1fK'%SampleList[-1]['Temperature'])
        if not G2frame.Contour:
            PlotList.reverse()
            ParmList.reverse()
            SampleList.reverse()
            LimitsList.reverse()
    if timeDebug:
        print('plot build time: %.3f for %dx%d patterns'%(time.time()-time0,len(PlotList[0][1][1]),len(PlotList)))
    lenX = 0  # length of first histogram, used for contour plots
    Ymax = None
    for ip,Pattern in enumerate(PlotList):
        xye = Pattern[1]
        xye = np.nan_to_num(xye)
        if xye[1] is None: continue
        if Ymax is None: Ymax = max(xye[1])
        Ymax = max(Ymax,max(xye[1]))
    if Ymax is None: return # nothing to plot
    offsetX = Page.plotStyle['Offset'][1]
    offsetY = Page.plotStyle['Offset'][0]
    if Page.plotStyle['logPlot']:
        Title = 'log('+Title+')'
    elif Page.plotStyle['sqrtPlot']:
        Title = r'$\sqrt{I}$ for '+Title
        Ymax = np.sqrt(Ymax)
    elif Page.plotStyle.get('WgtDiagnostic',False):
        Title = 'Scaling diagnostic for '+Title
    if G2frame.SubBack:
        Title += ' - background'
    if Page.plotStyle['qPlot'] or plottype in ['SASD','REFD'] and not G2frame.Contour and not ifLimits:
        xLabel = r'$Q, \AA^{-1}$'
    elif Page.plotStyle['dPlot'] and 'PWDR' in plottype and not ifLimits:
        xLabel = r'$d, \AA$'
    elif Page.plotStyle['chanPlot'] and G2frame.Contour:
        xLabel = 'Channel no.'
    else:
        if 'T' in ParmList[0]['Type'][0]:
            xLabel = r'$TOF, \mathsf{\mu}$s'
        elif 'E' in ParmList[0]['Type'][0]:
            xLabel = 'E, keV'
        else:
            xLabel = r'$\mathsf{2\theta}$'

    if G2frame.Weight and not G2frame.Contour:
        Plot.set_visible(False)         #hide old plot frame, will get replaced below
        GS_kw = {'height_ratios':[4, 1],}
        # try:
        Plot,Plot1 = Page.figure.subplots(2,1,sharex=True,gridspec_kw=GS_kw)
        # except AttributeError: # figure.Figure.subplots added in MPL 2.2
        #     Plot,Plot1 = MPLsubplots(Page.figure, 2, 1, sharex=True, gridspec_kw=GS_kw)
        Plot1.set_ylabel(r'$\mathsf{\Delta(I)/\sigma(I)}$',fontsize=16)
        Plot1.set_xlabel(xLabel,fontsize=16)
        Page.figure.subplots_adjust(left=16/100.,bottom=16/150.,
            right=.98,top=1.-16/200.,hspace=0)
    else:
        Plot.set_xlabel(xLabel,fontsize=16)
    if G2frame.Weight and G2frame.Contour:
        Title = r'$\mathsf{\Delta(I)/\sigma(I)}$ for '+Title
    if 'T' in ParmList[0]['Type'][0] or (Page.plotStyle['Normalize'] and not G2frame.SinglePlot):
        if Page.plotStyle['sqrtPlot']:
            Plot.set_ylabel(r'$\sqrt{Normalized\ intensity}$',fontsize=16)
        else:
            Plot.set_ylabel(r'$Normalized\ intensity$',fontsize=16)
    else:       #neutron TOF
        if 'PWDR' in plottype:
            if Page.plotStyle['sqrtPlot']:
                Plot.set_ylabel(r'$\sqrt{Intensity}$',fontsize=16)
            elif Page.plotStyle.get('WgtDiagnostic',False):
                Plot.set_ylabel('Intensity * Weight')
            elif G2frame.CumeChi and G2frame.SinglePlot:
                Plot.set_ylabel(r'$Intensity, cum('+Gkchisq+')$',fontsize=16)
            else:
                Plot.set_ylabel(r'$Intensity$',fontsize=16)
        elif plottype == 'SASD':
            if Page.plotStyle['sqPlot']:
                Plot.set_ylabel(r'$S(Q)=I*Q^{4}$',fontsize=16)
            else:
                Plot.set_ylabel(r'$Intensity,\ cm^{-1}$',fontsize=16)
        elif plottype == 'REFD':
            if Page.plotStyle['sqPlot']:
                Plot.set_ylabel(r'$S(Q)=R*Q^{4}$',fontsize=16)
            else:
                Plot.set_ylabel(r'$Reflectivity$',fontsize=16)                
    mpl.rcParams['image.cmap'] = G2frame.ContourColor
    mcolors = mpl.cm.ScalarMappable()       #wants only default as defined in previous line!!
    mcolors.set_array([]) # needed for MPL <=3.0.x
    if G2frame.Contour:
        ContourZ = []
        ContourY = []
        ContourX = None
        Xlist = []
        X0 = None
        Nseq = 0
    Nmax = len(PlotList)-1
    time0 = time.time()
    Plot.figure.subplots_adjust(right=.95)
    if G2frame.Contour and G2frame.TforYaxis:
        Plot.set_ylabel('Temperature',fontsize=14)
    elif G2frame.Contour:
        Plot.set_ylabel('Data sequence',fontsize=14)
    unequalArrays = False # set to True for contour plots with unequal pixels
    avgStep = None
    if G2frame.Contour:  # detect unequally spaced points in a contour plot
        for N,Pattern in enumerate(PlotList):
            xye = np.array(ma.getdata(Pattern[1])) # strips mask = X,Yo,W,Yc,Yb,Yd
            if Page.plotStyle['qPlot'] and 'PWDR' in plottype and not ifLimits:
                X = 2.*np.pi/G2lat.Pos2dsp(Parms,xye[0])
            elif Page.plotStyle['dPlot'] and 'PWDR' in plottype and not ifLimits:
                X = G2lat.Pos2dsp(Parms,xye[0])
            else:
                X = copy.deepcopy(xye[0])
            if not X0:   
                X0 = X[0] # save 1st point in 1st pattern
            elif abs(X0 - X[0]) > 0.05 * X0:
                unequalArrays = True
            if Page.plotStyle['qPlot'] or Page.plotStyle['dPlot']:  # not in original units
                unequalArrays = True
            elif 'T' in ParmList[0]['Type'][0] and not Page.plotStyle['chanPlot']: # assume TOF is non-linear steps
                unequalArrays = True
            # check to see if the average step size changes across the selected patterns
            elif avgStep is None and not unequalArrays:
                avgStep = (X[-1]-X[0])/(len(X)-1)
            elif not unequalArrays and abs(avgStep - (X[-1]-X[0])/(len(X)-1)) > 0.05 * avgStep:
                unequalArrays = True

    ExMask = []
    for N,Pattern in enumerate(PlotList):
        Parms = ParmList[N]
        Sample = SampleList[N]
        limits = np.array(LimitsList[N])
        ifpicked = False
        NoffY = offsetY*(Nmax-N)
        if Pattern[1] is None: continue # skip over uncomputed simulations
        xye = np.array(ma.getdata(Pattern[1])) # strips mask = X,Yo,W,Yc,Yb,Yd
        ExMask.append(np.full(len(xye[0]),False))
        if G2frame.PickId:   # when is this not true?
            ifpicked = Pattern[2] == G2frame.GPXtree.GetItemText(G2frame.PatternId)
            # recompute mask from excluded regions, in case they have changed
            xye0 = xye[0]  # no mask in case there are no limits
            for excl in limits[2:]:
                xye0 = ma.masked_inside(xye[0],excl[0],excl[1],copy=False)                   #excluded region mask
            if unequalArrays:
                xye0 = ma.masked_outside(xye[0],limits[1][0],limits[1][1],copy=False) #now mask for limits
                Lmask = ma.getmask(xye0)   # limits applied
                ExMask[N] = ExMask[N][~Lmask] # drop points outside limits
            elif not G2frame.Contour:
                xye0 = ma.masked_outside(xye0,limits[1][0],limits[1][1],copy=False) #now mask for limits
        else:
            xye0 = Pattern[1][0]  # keeps mask
            Lmask = Emask = np.full(len(xye0),False)

        if G2frame.Contour:
            xye0 = xye[0]   # drop mask for contouring

        # convert all X values and then reapply mask
        if Page.plotStyle['qPlot'] and 'PWDR' in plottype and not ifLimits:
            X = ma.array(2.*np.pi/G2lat.Pos2dsp(Parms,xye0.data),mask=xye0.mask)
        elif Page.plotStyle['dPlot'] and 'PWDR' in plottype and not ifLimits:
            X = ma.array(G2lat.Pos2dsp(Parms,xye0.data),mask=xye0.mask)
        else:
            X = copy.deepcopy(xye0)
        if ifpicked and not G2frame.Contour:
            savedX = copy.deepcopy(X)
            
        if not lenX:
            lenX = len(X)
        # show plot magnification factors
        magMarkers = []
        Plot.magLbls = []
        multArray = np.ones_like(Pattern[1][0])
        if 'PWDR' in plottype and G2frame.SinglePlot and not (
                Page.plotStyle['logPlot'] or Page.plotStyle['sqrtPlot'] or G2frame.Contour):
            if not refineMode:
                magLineList = data[0].get('Magnification',[])
            if ('C' in ParmList[0]['Type'][0] and Page.plotStyle['dPlot']) or (
                'B' in ParmList[0]['Type'][0] and Page.plotStyle['dPlot']) or (
                'E' in ParmList[0]['Type'][0] and Page.plotStyle['dPlot']) or (
                'T' in ParmList[0]['Type'][0] and Page.plotStyle['qPlot']): # reversed regions relative to data order
                tcorner = 1
                tpos = 1.0
                halign = 'right'
            else:
                tcorner = 0
                tpos = 1.0
                halign = 'left'
            ml0 = None
            for x,m in magLineList:
                ml = m
                if int(m) == m:
                    ml = int(m)
                if x is None:
                    magMarkers.append(None)
                    multArray *= m
                    ml0 = ml
                    continue
                multArray[Pattern[1][0]>x] = m
                if Page.plotStyle['qPlot']:
                    x = 2.*np.pi/G2lat.Pos2dsp(Parms,x)
                elif Page.plotStyle['dPlot']:
                    x = G2lat.Pos2dsp(Parms,x)
                # is range in displayed range (defined after newplot)?
                if not newPlot:
                    (xmin,xmax),ylim = G2frame.xylim
                    if x < xmin:
                        ml0 = ml
                        continue
                    if x > xmax:
                        continue
                # magnification region marker
                magMarkers.append(Plot.axvline(x,color='0.5',dashes=(1,1),
                                picker=True,pickradius=2.,label='_magline'))
                lbl = Plot.annotate("x{}".format(ml), xy=(x, tpos), xycoords=("data", "axes fraction"),
                    verticalalignment='bottom',horizontalalignment=halign,label='_maglbl')
                Plot.magLbls.append(lbl)
            if ml0:
                lbl = Plot.annotate("x{}".format(ml0), xy=(tcorner, tpos), xycoords="axes fraction",
                    verticalalignment='bottom',horizontalalignment=halign,label='_maglbl')
                Plot.magLbls.append(lbl)
                Page.toolbar.updateActions = (PlotPatterns,G2frame)
            multArray = ma.getdata(multArray)
        if 'PWDR' in plottype:
            YI = copy.copy(xye[1])      #yo
            if G2frame.SubBack:
                YI -= xye[4]            #background
            if Page.plotStyle['sqrtPlot']:
                olderr = np.seterr(invalid='ignore') #get around sqrt(-ve) error
                Y = np.where(YI>=0.,np.sqrt(YI),-np.sqrt(-YI))+NoffY*Ymax/100.0
                np.seterr(invalid=olderr['invalid'])
            elif Page.plotStyle.get('WgtDiagnostic',False):
                Y = xye[1]*xye[2]       #Y-obs*wt
            elif 'PWDR' in plottype and G2frame.SinglePlot and not (
                    Page.plotStyle['logPlot'] or Page.plotStyle['sqrtPlot'] or G2frame.Contour):
                Y = YI*multArray+NoffY*Ymax/100.0
            else:
                Y = YI+NoffY*Ymax/100.0
        elif plottype in ['SASD','REFD']:
            if plottype == 'SASD':
                B = xye[5]      #Yo - Yc
            else:
                B = np.zeros_like(xye[5])
            if Page.plotStyle['sqPlot']:
                Y = xye[1]*Sample['Scale'][0]*(1.05)**NoffY*X**4
            else:
                Y = xye[1]*Sample['Scale'][0]*(1.05)**NoffY
        if Page.plotStyle['exclude']:
            Y = ma.array(Y,mask=ma.getmask(X))
                
        if ifpicked and not G2frame.Contour: # draw limit & excluded region lines
            lims = limits[1:]
            if Page.plotStyle['qPlot'] and 'PWDR' in plottype and not ifLimits:
                lims = 2.*np.pi/G2lat.Pos2dsp(Parms,lims)
            elif Page.plotStyle['dPlot'] and 'PWDR' in plottype and not ifLimits:
                lims = G2lat.Pos2dsp(Parms,lims)
            # limit lines
            Lines.append(Plot.axvline(lims[0][0],color='g',dashes=(5,5),picker=True,pickradius=3.))    
            Lines.append(Plot.axvline(lims[0][1],color='r',dashes=(5,5),picker=True,pickradius=3.))
            # excluded region lines
            for i,item in enumerate(lims[1:]):
                Lines.append(Plot.axvline(item[0],color='m',dashes=(5,5),picker=True,pickradius=3.))    
                Lines.append(Plot.axvline(item[1],color='m',dashes=(5,5),picker=True,pickradius=3.))
                exclLines += [2*i+2,2*i+3]
        if G2frame.Contour:
            if Page.plotStyle['chanPlot']:
                if unequalArrays:
                    X = np.array(range(len(X)),float)
                else:
                    X = np.array(range(lenX),float)
                    Lmask = Emask = np.full(len(X),False)
            if G2frame.Weight:
                Ytmp = (xye[1]-xye[3])*np.sqrt(xye[2])
            else:
                Ytmp = Y
            # pad or truncate arrays when plotting with mpl.imshow
            if unequalArrays:
                ContourZ.append(ma.MaskedArray(Ytmp,Lmask).compressed())
            elif len(Y) < lenX:
                Yext = np.ones(lenX)*Ytmp[-1]
                Yext[:len(X)] = Ytmp
                ContourZ.append(Yext)
            elif len(Y) > lenX:
                ContourZ.append(Ytmp[:lenX])
            else:
                ContourZ.append(Ytmp)
            #if unequalArrays and G2frame.TforYaxis:
            #    TODO: could set this to temperature and then plot
            #    against temperature, but this only works if patterns are sorted by T
            ContourY.append(N)
            if unequalArrays:
                Xlist.append(ma.MaskedArray(X,Lmask).compressed())
            elif ContourX is None:
                ContourX = X
            Nseq += 1
        else:   # not contour plot
            if not G2frame.plusPlot:
                pP = ''
                lW = 1.5
            elif G2frame.plusPlot == 1:
                pP = '+'
                lW = 0
            else:
                pP = '+'
                lW = 1.5
            if plottype in ['SASD','REFD'] and Page.plotStyle['logPlot']:
                X *= (1.01)**(offsetX*N)
            else:
                xlim = Plot.get_xlim()
                DX = xlim[1]-xlim[0]
                X += 0.002*offsetX*DX*N
            if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Limits':
                Xum = ma.getdata(X) # unmasked version of X, use for data limits (only)
            else:
                Xum = X[:]
            if ifpicked:      # plotting of "master" pattern (histogram selected in data tree)
                ZI = copy.copy(xye[3])      #Yc
                if Page.plotStyle['sqrtPlot']:
                    olderr = np.seterr(invalid='ignore') #get around sqrt(-ve) error
                    Z = np.where(ZI>=0.,np.sqrt(ZI),-np.sqrt(-ZI))
                    np.seterr(invalid=olderr['invalid'])
                else:
                    if 'PWDR' in plottype and G2frame.SinglePlot and not (
                        Page.plotStyle['logPlot'] or Page.plotStyle['sqrtPlot'] or G2frame.Contour):
                        Z = ZI*multArray+NoffY*Ymax/100.0   #yc
                    else:
                        Z = ZI+NoffY*Ymax/100.0         #yc
                if 'PWDR' in plottype:
                    if Page.plotStyle['sqrtPlot']:
                        olderr = np.seterr(invalid='ignore') #get around sqrt(-ve) error
                        W = np.where(xye[4]>=0.,np.sqrt(xye[4]),-np.sqrt(-xye[4]))      #yb
                        np.seterr(invalid=olderr['invalid'])
                        D = np.where(xye[5],(Y-Z),0.)-Page.plotStyle['delOffset']
                    elif Page.plotStyle.get('WgtDiagnostic',False):
                        Z = D = W = np.zeros_like(Y)
                    elif G2frame.SinglePlot and not (
                        Page.plotStyle['logPlot'] or Page.plotStyle['sqrtPlot'] or G2frame.Contour):
                        W = xye[4]*multArray+NoffY*Ymax/100.0       #yb
                        D = multArray*xye[5]-Page.plotStyle['delOffset']  #Yo-Yc
                    else:
                        W = xye[4]+NoffY*Ymax/100.0     #yb
                        D = xye[5]-Page.plotStyle['delOffset']  #Yo-Yc
                elif plottype in ['SASD','REFD']:
                    if Page.plotStyle['sqPlot']:
                        W = xye[4]*X**4         #Yb
                        Z = xye[3]*X**4         #Yc
                        B = B*X**4              #(yo-Yc)*x**4
                    else:
                        W = xye[4]              #Yb
                    if G2frame.SubBack:
                        YB = Y-B
                        ZB = Z
                    else:
                        YB = Y
                        ZB = Z+B
                    try:
                        Plot.set_yscale("log",nonpositive='mask') # >=3.3
                    except:
                        Plot.set_yscale("log",nonpositive='mask')
                    if np.any(W>0.):
                        lims = [np.min(np.trim_zeros(W))/2.,np.max(Y)*2.]
                    else:
                        lims = [np.min(np.trim_zeros(YB))/2.,np.max(Y)*2.]
                    if all(np.isfinite(lims)): 
                        Plot.set_ylim(bottom=lims[0],top=lims[1])
                # Matplotlib artist lists used for refPlotUpdate
                ObsLine = None
                CalcLine = None
                BackLine = None
                DifLine = [None]
                if 'PWDR' in plottype and len(limits[2:]):   # compute mask for excluded regions
                    Emask = copy.deepcopy(ma.getmask(X))
                    for excl in limits[2:]:
                        Emask += ma.getmask(ma.masked_inside(xye[0],excl[0],excl[1],copy=False))
                    if Page.plotStyle['exclude']:            # optionally apply mask
                        Xum = ma.array(Xum,mask=Emask)
                        X = ma.array(X,mask=Emask)
                        Y = ma.array(Y,mask=Emask)
                        Z = ma.array(Z,mask=Emask)
                        W = ma.array(W,mask=Emask)
                    D = ma.array(D,mask=Emask)              # difference plot is always masked

                if G2frame.Weight:
                    Plot1.set_yscale("linear")                                                  
                    wtFactor = Pattern[0]['wtFactor']
                    if plottype in ['SASD','REFD']:
                        DZ = (Y-B-Z)*np.sqrt(wtFactor*xye[2])
                    else:
                        DZ = (xye[1]-xye[3])*np.sqrt(wtFactor*xye[2])
                        if 'PWDR' in plottype and len(limits[2:]):
                            DZ = ma.array(DZ,mask=Emask)   # weighted difference is always masked
                    DifLine = Plot1.plot(X,DZ,pwdrCol['Diff_color'],picker=True,pickradius=1.,label=incCptn('diff'))                    #(Io-Ic)/sig(Io)
                    Plot1.tick_params(labelsize=14)
                    Plot1.axhline(0.,color='k')
                    
                if Page.plotStyle['logPlot']:
                    if 'PWDR' in plottype:
                        try:
                            Plot.set_yscale("log",nonpositive='mask') # >=3.3
                        except:
                            Plot.set_yscale("log",nonpositive='mask')
                        Plot.plot(X,Y,marker=pP,color=pwdrCol['Obs_color'],linewidth=lW,picker=True,pickradius=3.,
                            clip_on=Clip_on,label=incCptn('obs'))
                        if G2frame.SinglePlot or G2frame.plusPlot:
                            Plot.plot(X,Z,pwdrCol['Calc_color'],picker=False,label=incCptn('calc'),linewidth=1.5)
                            if G2frame.plusPlot:
                                Plot.plot(X,W,pwdrCol['Bkg_color'],picker=False,label=incCptn('bkg'),linewidth=1.5)     #background
                    elif plottype in ['SASD','REFD']:
                        try:
                            Plot.set_xscale("log",nonpositive='mask') # >=3.3
                            Plot.set_yscale("log",nonpositive='mask')
                        except:
                            Plot.set_xscale("log",nonpositive='mask')
                            Plot.set_yscale("log",nonpositive='mask')
                        if G2frame.ErrorBars:
                            if Page.plotStyle['sqPlot']:
                                Plot.errorbar(X,YB,yerr=X**4*Sample['Scale'][0]*np.sqrt(1./(Pattern[0]['wtFactor']*xye[2])),
                                    ecolor=pwdrCol['Obs_color'],
                                picker=True,pickradius=3.,clip_on=Clip_on)
                            else:
                                Plot.errorbar(X,YB,yerr=Sample['Scale'][0]*np.sqrt(1./(Pattern[0]['wtFactor']*xye[2])),
                                    ecolor=pwdrCol['Obs_color'],
                                picker=True,pickradius=3.,clip_on=Clip_on,label=incCptn('obs'))
                        else:
                            Plot.plot(X,YB,marker=pP,color=pwdrCol['Obs_color'],linewidth=lW,
                                picker=True,pickradius=3.,clip_on=Clip_on,label=incCptn('obs'))
                        Plot.plot(X,W,pwdrCol['Calc_color'],picker=False,label=incCptn('bkg'),linewidth=1.5)     #const. background
                        Plot.plot(X,ZB,pwdrCol['Bkg_color'],picker=False,label=incCptn('calc'),linewidth=1.5)
                else:  # not logPlot
                    ymax = 1.
                    if Page.plotStyle['Normalize'] and Y.max() != 0 and not G2frame.SinglePlot:
                        ymax = Y.max()
                    if G2frame.SubBack:
                        if 'PWDR' in plottype:
                            ObsLine = Plot.plot(Xum,Y/ymax,color=pwdrCol['Obs_color'],marker=pP,linewidth=lW,
                                picker=False,clip_on=Clip_on,label=incCptn('obs-bkg'))  #Io-Ib
                            if np.any(Z):       #only if there is a calc pattern
                                CalcLine = Plot.plot(X,(Z-W)/ymax,pwdrCol['Calc_color'],picker=False,
                                    label=incCptn('calc-bkg'),linewidth=1.5)               #Ic-Ib
                        else:
                            Plot.plot(X,YB,color=pwdrCol['Obs_color'],marker=pP,linewidth=lW,
                                picker=True,pickradius=3.,clip_on=Clip_on,label=incCptn('obs'))
                            Plot.plot(X,ZB,pwdrCol['Bkg_color'],picker=False,label=incCptn('calc'),linewidth=1.5)
                    else:
                        if 'PWDR' in plottype:
                            ObsLine = Plot.plot(Xum,Y/ymax,color=pwdrCol['Obs_color'],marker=pP,linewidth=lW,
                                picker=True,pickradius=3.,clip_on=Clip_on,label=incCptn('obs'))    #Io
                            CalcLine = Plot.plot(X,Z/ymax,pwdrCol['Calc_color'],picker=False,label=incCptn('calc'),linewidth=1.5)                 #Ic
                        else:
                            Plot.plot(X,YB,color=pwdrCol['Obs_color'],marker=pP,linewidth=lW,
                                picker=True,pickradius=3.,clip_on=Clip_on,label=incCptn('obs'))
                            Plot.plot(X,ZB,pwdrCol['Bkg_color'],picker=False,label=incCptn('calc'),linewidth=1.5)
                    if 'PWDR' in plottype and (G2frame.SinglePlot and G2frame.plusPlot):
                        BackLine = Plot.plot(X,W/ymax,pwdrCol['Bkg_color'],picker=False,label=incCptn('bkg'),linewidth=1.5)                 #Ib
                        if not G2frame.Weight and np.any(Z):
                            DifLine = Plot.plot(X,D/ymax,pwdrCol['Diff_color'],linewidth=1.5,
                                picker=True,pickradius=1.,label=incCptn('diff'))                 #Io-Ic
                    Plot.axhline(0.,color='k',label='_zero')
                    
                    Plot.tick_params(labelsize=14)
                    # write a .csv file; not fully tested, but probably works where allowed
                    if 'PWDR' in plottype and G2frame.SinglePlot and plotOpt['saveCSV']:
                        plotOpt['saveCSV'] = False
                        fp = open(plotOpt['CSVfile'],'w')
                        G2plt.Write2csv(fp,['"limits"',lims[0][0],lims[0][1]])
                        l = []
                        PeakId = G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List')
                        peaks = G2frame.GPXtree.GetItemPyData(PeakId)
                        if inXtraPeakMode:
                            tbl = peaks['xtraPeaks']
                        else:
                            tbl = peaks['peaks']
                        for i,item in enumerate(tbl):
                            if type(item) is dict: continue
                            pos = item[0]
                            if Page.plotStyle['qPlot']:
                                l.append(2.*np.pi/G2lat.Pos2dsp(Parms,pos))
                            elif Page.plotStyle['dPlot']:
                                l.append(G2lat.Pos2dsp(Parms,pos))
                            else:
                                l.append(pos)
                        if l: G2plt.Write2csv(fp,['"peaks"']+l)
                        peaks['LaueFringe'] = peaks.get('LaueFringe',{})
                        l = []
                        for pos in peaks['LaueFringe'].get('satellites',[]):
                            if Page.plotStyle['qPlot']:
                                l.append(2.*np.pi/G2lat.Pos2dsp(Parms,pos))
                            elif Page.plotStyle['dPlot']:
                                l.append(G2lat.Pos2dsp(Parms,pos))
                            else:
                                l.append(pos)
                        if l: G2plt.Write2csv(fp,['"satellites"']+l)

                        G2plt.Write2csv(fp,['masked X','X','obs','calc','bkg','diff'],header=True)
                        for i in range(len(X)):
                            G2plt.Write2csv(fp,[X[i],X.data[i],Y[i],Z[i],W[i],D[i]],header=False)
                        fp.close()
                        print('file',plotOpt['CSVfile'],'written')
                        
                Page.SetToolTipString('')
                if G2frame.PickId:
                    if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Peak List':
                        tip = 'On data point: Pick peak - L or R MB. On line: L-move, R-delete'
                        Page.SetToolTipString(tip)
                        peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Peak List'))
                        if inXtraPeakMode:
                            peaks['xtraPeaks'] = peaks.get('xtraPeaks',[])
                            tbl = peaks['xtraPeaks']
                            color = 'r'
                        else:
                            tbl = peaks['peaks']
                            color = 'b'
                        try:
                            # find any peak rows that are selected
                            selectedPeaks = list(set(
                                [row for row,col in G2frame.reflGrid.GetSelectedCells()] +
                                G2frame.reflGrid.GetSelectedRows()))
#                            G2frame.dataWindow.movePeak.Enable(len(selectedPeaks) == 1) # allow peak move from table when one peak is selected
                            for i,item in enumerate(tbl):
                                if type(item) is dict: continue
                                if i in selectedPeaks:
                                    Ni = N+1
                                    plotVline(Page,Plot,Lines,Parms,item[0],'yellow',False,'-')
                                    Lines[-1].set_lw(Lines[-1].get_lw()+1)
                                    plotVline(Page,Plot,Lines,Parms,item[0],color,True)
                                    Lines[-1].set_lw(Lines[-1].get_lw()+1)
                                else:
                                    Ni = N
                                    plotVline(Page,Plot,Lines,Parms,item[0],color,True)
                        except:
                            pass
                        peaks['LaueFringe'] = peaks.get('LaueFringe',{})
                        SatLines = []
                        for pos in peaks['LaueFringe'].get('satellites',[]):
                            plotVline(Page,Plot,SatLines,Parms,pos,'k',False)
#                        for pos in peaks['xtraPeaks']:
#                            plotVline(Page,Plot,Lines,Parms,pos[0],'r',False)
                    if G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Limits':
                        tip = 'On data point: Lower limit - L MB; Upper limit - R MB. On limit: MB down to move'
                        Page.SetToolTipString(tip)
                        limits = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Limits'))  # used anywhere?
                        
            else:   #not picked
                ymax = 1.
                if Page.plotStyle['Normalize'] and Y.max() != 0:
                    ymax = Y.max()
                icolor = 256*N//len(PlotList)
                if Page.plotStyle['logPlot']:
                    if 'PWDR' in plottype:
                        try:
                            Plot.semilogy(X,Y,color=mcolors.cmap(icolor), # >=3.3
                                picker=False,nonpositive='mask',linewidth=1.5)
                        except:
                            Plot.semilogy(X,Y,color=mcolors.cmap(icolor),
                                picker=False,nonpositive='mask')
                    elif plottype in ['SASD','REFD']:
                        try:
                            Plot.semilogy(X,Y,color=mcolors.cmap(icolor),
                                picker=False,nonpositive='mask',linewidth=1.5)
                        except:
                            Plot.semilogy(X,Y,color=mcolors.cmap(icolor),
                                picker=False,nonpositive='mask')
                else:
                    if 'PWDR' in plottype:
                        Plot.plot(X,Y/ymax,color=mcolors.cmap(icolor),picker=False)
                    elif plottype in ['SASD','REFD']:
                        try:
                            Plot.loglog(X,Y,mcolors.cmap(icolor),
                                picker=False,nonpositive='mask',linewidth=1.5)
                        except:
                            Plot.loglog(X,Y,mcolors.cmap(icolor),
                                picker=False,nonpositive='mask')
                        Plot.set_ylim(bottom=np.min(np.trim_zeros(Y))/2.,top=np.max(Y)*2.)
                            
                if Page.plotStyle['logPlot'] and 'PWDR' in plottype:
                    Plot.set_ylim(bottom=np.min(np.trim_zeros(Y))/2.,top=np.max(Y)*2.)
    #============================================================
    # plot HKL labels
    data[0]['HKLconfig'] = data[0].get('HKLconfig',{})
    alpha = data[0]['HKLconfig'].get('alpha',2)
    font = int(data[0]['HKLconfig'].get('Font','8'))
    angle = int(90 * data[0]['HKLconfig'].get('Orientation',0))
    props = dict(boxstyle='round,pad=0.15', facecolor='#ccc',
                     alpha=(4 - alpha)/4., ec='#ccc')
    markHKLs = (G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Reflection Lists' or 
        'PWDR' in G2frame.GPXtree.GetItemText(G2frame.PickId) or refineMode
        or (inXtraPeakMode and G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Peak List'))
    for ph,markerlist in data[0].get('HKLmarkers',{}).items():
        if Page.plotStyle['logPlot'] or not markHKLs: continue
        if ph not in Page.phaseColors: continue
        color = Page.phaseColors[ph]
        for key,(ypos,hkls) in markerlist.items():
            if len(hkls) == 0: continue
            comma = False             # are commas needed?
            for hkl in hkls:
                if np.any([True for i in hkl if abs(i) > 9]):
                    comma = True
                    break
            lbl = ''
            for hkl in hkls:
                hkllbl = ''
                if lbl: lbl += '\n'
                for i in hkl:
                    if hkllbl and comma: hkllbl += ','
                    if i < 0:
                        lbl += r'$\overline{' + f'{-i}' + r'}$'
                    else:
                        lbl += f'{i}'
            xpos = float(key)
            if Page.plotStyle['qPlot']:
                xpos = 2.*np.pi/G2lat.Pos2dsp(Parms,xpos)
            elif Page.plotStyle['dPlot']:
                xpos = G2lat.Pos2dsp(Parms,xpos)
            if Page.plotStyle['sqrtPlot']:
                ypos = np.sqrt(abs(ypos))*np.sign(ypos)
            artist = Plot.text(xpos,ypos,lbl,fontsize=font,c=color,ha='center',
                      va='top',bbox=props,picker=True,rotation=angle,
                      label='_'+ph)
            artist.key = (ph,key)
    #============================================================
    if timeDebug:
        print('plot fill time: %.3f'%(time.time()-time0))
    if Page.plotStyle.get('title',True) and not magLineList:
        Plot.set_title(Title)

    if G2frame.PickId and not G2frame.Contour:
        Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
        orange = [255/256.,128/256.,0.]
        if 'PWDR' in plottype and G2frame.SinglePlot and G2frame.CumeChi:
            CY = np.cumsum(W**2*(Y-Z)**2)
            scale = np.max(CY)/np.max(Y)
            CY /= scale
            Plot.plot(X,CY,'k',picker=False,label='cum('+Gkchisq+')')
        selectedPeaks = []
        if G2frame.GPXtree.GetItemText(G2frame.PickId) in ['Index Peak List']:
            # find any peak rows that are selected
            selectedPeaks = list(set(
                [row for row,col in G2frame.indxPeaks.GetSelectedCells()] +
                    G2frame.indxPeaks.GetSelectedRows()))
        if G2frame.GPXtree.GetItemText(G2frame.PickId) in ['Index Peak List','Unit Cells List']:
            peaks = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Index Peak List'))
            if len(peaks):  # are there any peaks?
                for i,peak in enumerate(peaks[0]):
                    if Page.plotStyle['qPlot']:
                        x = 2.*np.pi/G2lat.Pos2dsp(Parms,peak[0])
                    elif Page.plotStyle['dPlot']:
                        x = G2lat.Pos2dsp(Parms,peak[0])
                    else:
                        x = peak[0]
                    if i in selectedPeaks:
                        Plot.axvline(x,color='yellow',lw=2)
                    if peak[2]:
                        if i in selectedPeaks:
                            Plot.axvline(x,color='yellow',lw=2)
                            Plot.axvline(x,color='b',lw=2,ls='dotted')
                        else:
                            Plot.axvline(x,color='b')
            for hkl in G2frame.HKL:
                clr = orange
                dash = (3,3)
                if len(hkl) > 6 and hkl[3]:
                    clr = 'g'
                hklind = G2frame.PlotOpts.get('hklHighlight',0)
                if hklind != 0:  # highlight selected classes of reflections
                    if hkl[hklind-1] != 0:
                        clr = 'b'
                        dash = (5,2)
                if Page.plotStyle['qPlot']:
                    Plot.axvline(2.*np.pi/G2lat.Pos2dsp(Parms,hkl[-2]),color=clr,dashes=dash,lw=1.5)
                elif Page.plotStyle['dPlot']:
                    Plot.axvline(G2lat.Pos2dsp(Parms,hkl[-2]),color=clr,dashes=dash,lw=1.5)
                else:
                    Plot.axvline(hkl[-2],color=clr,dashes=dash,lw=1.5)
            for hkl in G2frame.Extinct: # plot extinct reflections
                clr = 'b'
                if Page.plotStyle['qPlot']:
                    Plot.axvline(2.*np.pi/G2lat.Pos2dsp(Parms,hkl[-2]),color=clr,dashes=(3,3),lw=2)
                elif Page.plotStyle['dPlot']:
                    Plot.axvline(G2lat.Pos2dsp(Parms,hkl[-2]),color=clr,dashes=(3,3),lw=2)
                else:
                    Plot.axvline(hkl[-2],color=clr,dashes=(3,3),lw=2)
        elif Page.plotStyle.get('WgtDiagnostic',False):
            pass # skip reflection markers
        elif (G2frame.GPXtree.GetItemText(G2frame.PickId) in ['Reflection Lists','Limits']
                  or 'PWDR' in G2frame.GPXtree.GetItemText(G2frame.PickId)
                  or refineMode
                  or (inXtraPeakMode and
                    G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Peak List')
                ):
            l = GSASIIpath.GetConfigValue('Tick_length',8.0)
            w = GSASIIpath.GetConfigValue('Tick_width',1.)
            for pId,phase in enumerate(Page.phaseList):
                if 'list' in str(type(Phases[phase])):
                    continue
                if phase in Page.phaseColors:
                    plcolor = Page.phaseColors[phase]
                else: # how could this happen? 
                    plcolor = 'k'
                    #continue
                peaks = Phases[phase].get('RefList',[])
                if not len(peaks):
                    continue
                if Phases[phase].get('Super',False):
                    peak = np.array([[peak[5],peak[6]] for peak in peaks])
                else:
                    peak = np.array([[peak[4],peak[5]] for peak in peaks])
                pos = Page.plotStyle['refOffset']-pId*Page.plotStyle['refDelt']*np.ones_like(peak)
                if Page.plotStyle['qPlot']:
                    xtick = 2*np.pi/peak.T[0]
                elif Page.plotStyle['dPlot']:
                    xtick = peak.T[0]
                else:
                    xtick = peak.T[1]
                if not Page.plotStyle.get('flTicks',False):     # short tick-marks
                    Page.tickDict[phase],_ = Plot.plot(
                        xtick,pos,'|',mew=w,ms=l,picker=True,pickradius=3.,
                        label=phase,color=plcolor)
                    # N.B. above creates two Line2D objects, 2nd is ignored.
                    # Not sure what each does.
                else:                                           # full length tick-marks
                    if len(xtick) > 0:
                        # create an ~hidden tickmark to create a legend entry
                        Page.tickDict[phase] = Plot.plot(xtick[0],0,'|',mew=0.5,ms=l,
                                                label=phase,color=plcolor)[0]
                    for xt in xtick: # a separate line for each reflection position
                            Plot.axvline(xt,color=plcolor,
                                        picker=True,pickradius=3.,
                                        label='_FLT_'+phase,lw=0.5)
            handles,legends = Plot.get_legend_handles_labels() 
            if handles:
                labels = dict(zip(legends,handles))     # remove duplicate phase entries
                handles = [labels[item] for item in labels]
                legends = list(labels.keys())
                if len(Phases) and obsInCaption: 
                    Plot.legend(handles,legends,title='Phases & Data',loc='best')
                else:
                    Plot.legend(handles,legends,title='Data',loc='best')
    
    if G2frame.Contour:
        time0 = time.time()
        acolor = G2plt.GetColorMap(G2frame.ContourColor)
        Vmin = Ymax*G2frame.Cmin
        Vmax = Ymax*G2frame.Cmax
        if unequalArrays:
            if G2frame.Weight:  
                #Vmin = min([i.min() for i in ContourZ])
                Vmin = min([ma.array(i,mask=m).min() for i,m in zip(ContourZ,ExMask)]) # don't count excluded points in limits
                #Vmax = max([i.max() for i in ContourZ])
                Vmax = max([ma.array(i,mask=m).max() for i,m in zip(ContourZ,ExMask)])
            if G2frame.TforYaxis:
                imgLbls = Temps
            else:
                imgLbls = []
            uneqImgShow(Plot.figure,Plot,Xlist,ContourZ,cmap=acolor,
                         vmin=Vmin,vmax=Vmax,Ylbls=imgLbls)
            Page.Img = None   # don't have an overall image
            if G2frame.TforYaxis:
                Plot.yaxis.set_label_coords(-.1, .5)
            else:
                Plot.yaxis.set_label_coords(-.05, .5)
            Plot.xaxis.set_label_coords(0.5, -.07)
        else:
            if G2frame.Weight:
                Vmin = np.min(ContourZ)
                Vmax = np.max(ContourZ)
            Page.Img = Plot.imshow(ContourZ,cmap=acolor,vmin=Vmin,vmax=Vmax,
                interpolation=G2frame.Interpolate,extent=[ContourX[0],ContourX[-1],ContourY[0]-.5,ContourY[-1]+.5],
                aspect='auto',origin='lower')
            if G2frame.TforYaxis:
                imgAx = Page.Img.axes
                ytics = imgAx.get_yticks()
                # ytics = np.where(ytics<len(Temps),ytics,-1)
                # imgAx.set_yticks(ytics)
                ylabs = [Temps[int(i)] for i in ytics[:-1]]
                imgAx.set_yticklabels(ylabs)
            Page.figure.colorbar(Page.Img)
        if timeDebug:
            print('Contour display time: %.3f'%(time.time()-time0))
    else:
        G2frame.Lines = Lines
        G2frame.MagLines = magMarkers
    if G2frame.PickId and G2frame.GPXtree.GetItemText(G2frame.PickId) == 'Background':
        mag2th = [0]+[x for x,m in data[0].get('Magnification',[])][1:]
        magmult = [m for x,m in data[0].get('Magnification',[])]
        # plot fixed background points
        backDict = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Background'))[1]
        try:
            Parms,Parms2 = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PatternId, 'Instrument Parameters'))
        except TypeError:
            Parms = None
        for x,y in backDict.get('FixedPoints',[]):
            if magmult:
                mult = magmult[np.searchsorted(mag2th, x, side = 'right')-1]
            else:
                mult = 1.
            if G2frame.Weight:
                axis = Page.figure.axes[1]
            else:
                axis = Page.figure.gca()
            # "normal" intensity modes only!
            if G2frame.SubBack or G2frame.Weight or G2frame.Contour or not G2frame.SinglePlot:
                break
            if y < 0 and (Page.plotStyle['sqrtPlot'] or Page.plotStyle['logPlot']):
                y = axis.get_ylim()[0] # put out of range point at bottom of plot
            elif Page.plotStyle['sqrtPlot']:
                y = math.sqrt(y)
            if Page.plotStyle['qPlot']:     #Q - convert from 2-theta
                if Parms:
                    x = 2*np.pi/G2lat.Pos2dsp(Parms,x)
                else:
                    break
            elif Page.plotStyle['dPlot']:   #d - convert from 2-theta
                if Parms:
                    x = G2lat.Dsp2pos(Parms,x)
                else:
                    break
            Plot.plot(x,y*mult,'rD',clip_on=Clip_on,picker=True,pickradius=10.)

    # plot the partials
    plotOpt['lineList']  = ['obs','calc','bkg','zero','diff']
    if 'PWDR' in plottype and G2frame.SinglePlot and Page.plotStyle['partials'] and 'hId' in data[0]:
        initPartialOpts(Page.phaseColors)
        x, yb, ypList = G2frame.LoadPartial(data[0]['hId'])            
        if x is not None and len(ypList) > 1:
            if Page.plotStyle['qPlot']:
                x = 2.*np.pi/G2lat.Pos2dsp(Parms,x)
            elif Page.plotStyle['dPlot']:
                x = G2lat.Pos2dsp(Parms,x)
            olderr = np.seterr(invalid='ignore') #get around sqrt(-ve) error
            for ph in ypList:
                if not partialOpts[ph]['Show']: continue
                pcolor = partialOpts[ph]['color']
                pwidth = partialOpts[ph]['width']
                pLinStyl = partialOpts[ph]['MPLstyle']
                if G2frame.SubBack:
                    y = ypList[ph]
                else:
                    y = ypList[ph]+yb
                if Page.plotStyle['sqrtPlot']:
                    y = np.where(y>=0.,np.sqrt(y),-np.sqrt(-y))
                Plot.plot(x,y,pcolor,picker=False,label=ph,linewidth=pwidth,
                              linestyle=pLinStyl)
    if not newPlot:
        # this restores previous plot limits (but I'm not sure why there are two .push_current calls)
        Page.toolbar.push_current()
        if G2frame.Contour: # for contour plots expand y-axis to include all histograms
            G2frame.xylim = (G2frame.xylim[0], (0.,len(PlotList)))
        if 'PWDR' in plottype:
            Plot.set_xlim(G2frame.xylim[0])
            Plot.set_ylim(G2frame.xylim[1])
        Page.toolbar.push_current()
        Page.ToolBarDraw()
    else:
        G2frame.xylim = Plot.get_xlim(),Plot.get_ylim()
        Page.canvas.draw()
    olderr = np.seterr(invalid='ignore') #ugh - this removes a matplotlib error for mouse clicks in log plots
    # and sqrt(-ve) in np.where usage               
    if 'PWDR' in G2frame.GPXtree.GetItemText(G2frame.PickId):
        if len(Page.tickDict.keys()) == 1:
            G2frame.dataWindow.moveTickLoc.Enable(True)
        elif len(Page.tickDict.keys()) > 1:
            G2frame.dataWindow.moveTickLoc.Enable(True)
            G2frame.dataWindow.moveTickSpc.Enable(True)
        if DifLine[0]:
            G2frame.dataWindow.moveDiffCurve.Enable(True)
    if refineMode: return refPlotUpdate

def PublishRietveldPlot(G2frame,Pattern,Plot,Page,reuse=None):
    '''Creates a window to show a customizable "Rietveld" plot. Exports that 
    plot as a publication-quality file. Will only work only when a single 
    pattern is displayed.

    :param wx.Frame G2Frame: the main GSAS-II window
    :param list Pattern: list of np.array items with obs, calc (etc.) diffraction pattern
    :param mpl.axes Plot: axes of the graph in plot window
    :param wx.Panel Page: tabbed panel containing the plot
    '''
    def Init_fmts():
        figure = mplfig.Figure(dpi=200,figsize=(6,8))
        canvas = hcCanvas(figure)
        fmtDict = canvas.get_supported_filetypes()
        plotWidgets['fmtChoices'] = [fmtDict[j]+', '+j for j in sorted(fmtDict)]
        plotWidgets['fmtChoices'].append('Data file with plot elements, csv')
        plotWidgets['fmtChoices'].append('Grace input file, agr')
        plotWidgets['fmtChoices'].append('Igor Pro input file, itx')
        if sys.platform == "win32":
            plotWidgets['fmtChoices'].append('OriginPro connection')
    def Initialize():
        '''Set up initial values in plotOpt
        '''
        plotOpt['initNeeded'] = False
        # create a temporary hard-copy figure to get output options
        figure = mplfig.Figure(dpi=200,figsize=(6,8))
        canvas = hcCanvas(figure)
        figure.clear()
        fmtDict = canvas.get_supported_filetypes()
        if plotOpt['format'] is None:
            if 'pdf' in fmtDict:
                plotOpt['format'] = fmtDict['pdf'] + ', pdf'
            else:
                plotOpt['format'] = plotWidgets['fmtChoices'][0]
        plotOpt['lineWid'] = '1'
        plotOpt['tickSiz'] = '6'
        plotOpt['tickWid'] = '1'
        plotOpt['markerWid'] = '1'
        plotOpt['markerSiz'] = '8'
        plotOpt['markerSym'] = '+'
        lims = {}
        lims['xmin'],lims['xmax'] = Plot.get_xlim()
        lims['ymin'],lims['ymax'] = Plot.get_ylim()
        for key in 'xmin','xmax','ymin','ymax':
            if not plotOpt[key+'_use']: plotOpt[key] = lims[key]
        #if mpl.__version__.split('.')[0] == '1':
        #    G2G.G2MessageBox(G2frame.plotFrame,
        #        ('You are using an older version of Matplotlib ({}). '.format(mpl.__version__) +
        #        '\nPlot quality will be improved by updating (use conda update matplotlib)'),
        #                     'Old matplotlib')
        
    def GetColors():
        '''Set up initial values in plotOpt for colors and legend. Also catalog phases.
        '''
        if hasattr(mpcls,'to_rgba'):
            MPL2rgba = mpcls.to_rgba
        else:
            MPL2rgba = mpcls.ColorConverter().to_rgba
        for i,l in enumerate(Plot.lines):
            lbl = l.get_label()
            if 'magline' in lbl:
                pass
            elif lbl[1:] in plotOpt['lineList']: # item not in legend
                if lbl[1:] in plotOpt['colors']: continue
                plotOpt['colors'][lbl[1:]] = MPL2rgba(l.get_color())
                plotOpt['legend'][lbl[1:]] = False
            elif lbl in plotOpt['lineList']:
                if lbl in plotOpt['colors']: continue
                plotOpt['colors'][lbl] = MPL2rgba(l.get_color())
                plotOpt['legend'][lbl] = True
        plotWidgets['phaseList'].clear()
        for ph,l in Page.tickDict.items():
            plotWidgets['phaseList'].append(ph)
            if lbl in plotOpt['colors']: continue
            plotOpt['colors'][ph] = MPL2rgba(l.get_color())
            plotOpt['legend'][ph] = True
        return

    def RefreshPlot(*args,**kwargs):
        '''Update the plot on the dialog
        '''
        figure.clear()
        CopyRietveldPlot(G2frame,Pattern,Plot,Page,figure,plotWidgets['phaseList'])
        figure.canvas.draw()
        
        
    def CopyRietveld2Grace(Pattern,Plot,Page,plotOpt,filename):
        '''Copy the contents of the Rietveld graph from the plot window to
        a Grace input file (tested with QtGrace). Uses values from Pattern 
        to also generate a delta/sigma plot below. 
        '''
        def ClosestColorNumber(color):
            '''Convert a RGB value to the closest default Grace color
            '''
            import matplotlib.colors as mpcls
            colorlist = ('white','black','red','green','blue','yellow','brown',
                            'gray','purple','cyan','magenta','orange') # ordered by grace's #
            if not hasattr(mpcls,'to_rgba'): mpcls = mpcls.ColorConverter()
            return (np.sum(([np.array(mpcls.to_rgb(c)) for c in colorlist] -
                                np.array(color[:3]))**2,axis=1)).argmin()

        # blocks of code used in grace .agr files
        linedef = '''@{0} legend "{1}"
@{0} line color {2}
@{0} errorbar color {2}
@{0} symbol color {2}
@{0} symbol {3}
@{0} symbol fill color {2}
@{0} linewidth {4}
@{0} symbol linewidth {6}
@{0} line type {7}
@{0} symbol size {5}
@{0} symbol char 46
@{0} symbol fill pattern 1
@{0} hidden false
@{0} errorbar off\n'''
        linedef1 = '''@{0} legend "{1}"
@{0} line color {2}
@{0} errorbar color {2}
@{0} symbol color {2}
@{0} symbol fill color {2}
@{0} symbol 11
@{0} linewidth 0
@{0} linestyle 0
@{0} symbol size {3}
@{0} symbol linewidth {4}
@{0} symbol char 124
@{0} symbol fill pattern 1
@{0} hidden false\n'''
        linedef2 = '''@{0} legend "{1}"
@{0} line color {2}
@{0} errorbar color {2}
@{0} symbol color {2}
@{0} symbol fill color {2}
@{0} symbol 0
@{0} linewidth 0
@{0} linestyle 0
@{0} symbol size 1
@{0} symbol linewidth 0
@{0} symbol char 124
@{0} symbol fill pattern 1
@{0} hidden false
@{0} errorbar on
@{0} errorbar size 0
@{0} errorbar riser linewidth {3}\n'''
        linedef3 = '''@{0} legend "{1}"
@{0} line color {2}
@{0} errorbar color {2}
@{0} symbol color {2}
@{0} symbol {3}
@{0} symbol fill color {2}
@{0} linewidth {4}
@{0} symbol linewidth {6}
@{0} line type {7}
@{0} symbol size {5}
@{0} symbol char 46
@{0} symbol fill pattern 1
@{0} hidden false
@{0} errorbar off\n'''    
        grace_symbols = {"":0, "o":1 ,"s":2, "D":3, "^":4, "3":5, 'v':6,
            "4": 7, "+":8, "P":8, "x":9, "X":9, "*":10, ".":11}

        fp = open(filename,'w')
        fp.write("# Grace project file\n#\n@version 50010\n")
        # size of plots on page
        xmar = (.15,1.2)
        ymar = (.15,.9)
        top2bottom = 4. # 4 to 1 spacing for top to bottom boxes

        # scaling for top box
        fp.write('@g{0} hidden false\n@with g{0}\n@legend {1}\n'.format(0,"on"))
        fp.write('@legend {}, {}\n'.format(xmar[1]-.2,ymar[1]-.05))
        fp.write('@world xmin {}\n@world xmax {}\n'.format(Plot.get_xlim()[0],Plot.get_xlim()[1]))
        fp.write('@world ymin {}\n@world ymax {}\n'.format(Plot.get_ylim()[0],Plot.get_ylim()[1]))
        fp.write('@view xmin {}\n@view xmax {}\n'.format(xmar[0],xmar[1]))
        fp.write('@view ymin {}\n@view ymax {}\n'.format((1./top2bottom)*(ymar[1]-ymar[0])+ymar[0],ymar[1]))
        xticks = Plot.get_xaxis().get_majorticklocs()
        fp.write('@{}axis tick major {}\n'.format('x',xticks[1]-xticks[0]))
        yticks = Plot.get_yaxis().get_majorticklocs()    
        fp.write('@{}axis tick major {}\n'.format('y',yticks[1]-yticks[0]))
        fp.write('@{}axis ticklabel char size {}\n'.format('x',0)) # turns off axis labels
        if 'sqrt' in Plot.yaxis.get_label().get_text():
            #ylbl = 'sqrt(Intensity)' # perhaps there is a way to get the symbol in xmgrace but I did not find it
            ylbl = r'\x\#{d6}\f{}\oIntensity\O' # from Carlo Segre
        else:
            ylbl = 'Intensity'
        fp.write('@{0}axis label "{1}"\n@{0}axis label char size {2}\n'.format(
            'y',ylbl,float(plotOpt['labelSize'])/8.))
        fp.write('@{0}axis label place spec\n@{0}axis label place {1}, {2}\n'.format('y',0.0,0.1))
    # ======================================================================
    # plot magnification lines and labels (first, so "under" data)
        for i,l in enumerate(Plot.lines):
            lbl = l.get_label()
            if 'magline' not in lbl: continue
            #ax0.axvline(l.get_data()[0][0],color='0.5',dashes=(1,1))
            # vertical line
            s = '@with line\n@ line on\n@ line loctype world\n@ line g0\n'
            fp.write(s)
            s = '@ line {0}, {1}, {0}, {2}\n'
            fp.write(s.format(
                l.get_data()[0][0],Plot.get_ylim()[0],Plot.get_ylim()[1]))
            s = '@ line linewidth 2\n@ line linestyle 2\n@ line color 1\n@ line arrow 0\n@line def\n'
            fp.write(s)
        for l in Plot.texts:
            if 'magline' not in l.get_label(): continue
            if l.xycoords[0] == 'data':
                xpos = l.get_position()[0]
            elif l.get_position()[0] == 0:
                xpos = Plot.get_xlim()[0]
            else:
                xpos = Plot.get_xlim()[1]*.95
            s = '@with string\n@    string on\n@    string loctype world\n'
            fp.write(s)
            s = '@    string g{0}\n@    string {1}, {2}\n@    string color 1\n'
            fp.write(s.format(0,xpos,Plot.get_ylim()[1]))
            s = '@    string rot 0\n@    string font 0\n@    string just 4\n'
            fp.write(s)
            s = '@    string char size {1}\n@    string def "{0}"\n'
            fp.write(s.format(l.get_text(),float(plotOpt['labelSize'])/8.))
        datnum = -1
    # ======================================================================
    # plot data 
        for l in Plot.lines:
            if l.get_label() in ('obs','calc','bkg','zero','diff'):
                lbl = l.get_label()
            elif l.get_label()[1:] in ('obs','calc','bkg','zero','diff'):
                lbl = l.get_label()[1:]
            else:
                continue
            c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
            gc = ClosestColorNumber(c)
            if sum(c) == 4.0: # white on white, skip
                continue
            marker = l.get_marker()
            lineWid = l.get_lw()
            siz = l.get_markersize()
            mkwid = l.get_mew()
            glinetyp = 1
            if lbl == 'obs':
                gsiz = float(plotOpt['markerSiz'])/8.
                marker = plotOpt['markerSym']
                gmw = float(plotOpt['markerWid'])
                gsym = grace_symbols.get(marker,5)
                glinetyp = 0
            else:
                gsym = 0
                gsiz = 0
                gmw = 0
                lineWid = float(plotOpt['lineWid'])
            if not plotOpt['Show'].get(lbl,True): continue
            if plotOpt['legend'].get(lbl):
                glbl = lbl
            else:
                glbl = ""
            datnum += 1
            fp.write("@type xy\n")
            if lbl == 'zero':
                fp.write("{} {}\n".format(Plot.get_xlim()[0],0))
                fp.write("{} {}\n".format(Plot.get_xlim()[1],0))
            elif not ma.any(l.get_xdata().mask):
                for x,y in zip(l.get_xdata(),l.get_ydata()):
                    fp.write("{} {}\n".format(x,y))
            else:
                for x,y,m in zip(l.get_xdata(),l.get_ydata(),l.get_xdata().mask):
                    if not m: fp.write("{} {}\n".format(x,y))
            fp.write("&\n")
            fp.write(linedef.format("s"+str(datnum),glbl,gc,gsym,lineWid,gsiz,gmw,glinetyp))
    #======================================================================
    # reflection markers. Create a single hidden entry for the legend
    # and use error bars for the vertical lines
        for l in Plot.lines:
            glbl = lbl = l.get_label()
            if l not in Page.tickDict.values(): continue
            c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
            gc = ClosestColorNumber(c)
            siz = float(plotOpt['tickSiz'])*(Plot.get_ylim()[1] - Plot.get_ylim()[0])/(100*6) # 1% for siz=6
            mkwid = float(plotOpt['tickWid'])
            if sum(c) == 4.0: continue # white: ignore
            if not plotOpt['Show'].get(lbl,True): continue
            if plotOpt['legend'].get(lbl):
                # invisible data point for 
                datnum += 1
                fp.write("@type xy\n")
#                fp.write("-1 -1\n".format(x,y))
                fp.write("-1 -1\n")
                fp.write(linedef1.format(
                    "s"+str(datnum),glbl,gc,float(plotOpt['tickSiz'])/8.,mkwid))
            # plot values with error bars
            datnum += 1
            fp.write("@type xydy\n")
            for x,y in zip(l.get_xdata(),l.get_ydata()):
                fp.write("{} {} {}\n".format(x,y,siz))
            fp.write("&\n")
            fp.write(linedef2.format("s"+str(datnum),'',gc,mkwid))
    #======================================================================
    # Start (obs-cal)/sigma plot
        wtFactor = Pattern[0]['wtFactor']
        ysig = (Pattern[1][1]-Pattern[1][3])*np.sqrt(wtFactor*Pattern[1][2])
        fp.write("@type xy\n")
        # scaling for bottom box
        fp.write('@g{0} hidden false\n@with g{0}\n@legend {1}\n'.format(1,"off"))
        fp.write('@world xmin {}\n@world xmax {}\n'.format(Plot.get_xlim()[0],Plot.get_xlim()[1]))
        fp.write('@world ymin {}\n@world ymax {}\n'.format(ysig.min(),ysig.max()))
        fp.write('@view xmin {}\n@view xmax {}\n'.format(xmar[0],xmar[1]))
        fp.write('@view ymin {}\n@view ymax {}\n'.format(
            ymar[0],(1./top2bottom)*(ymar[1]-ymar[0])+ymar[0]))
        if 'theta' in Plot.get_xlabel():
            xlbl = r'2\f{Symbol}q'
        elif 'TOF' in Plot.get_xlabel():
            xlbl = r'TOF, \f{Symbol}m\f{}s'
        else:
            xlbl = Plot.get_xlabel().replace('$','')        
        fp.write('@{0}axis label "{1}"\n@{0}axis label char size {2}\n'.format(
            'x',xlbl,float(plotOpt['labelSize'])/8.))
        fp.write('@{0}axis label "{1}"\n@{0}axis label char size {2}\n'.format(
            'y',r'\f{Symbol}D/s',float(plotOpt['labelSize'])/8.))
        xticks = Plot.get_xaxis().get_majorticklocs()
        # come up with a "nice" tick interval for (o-c)/sig, since I am not sure
        # if this can be defaulted
        ytick = (ysig.max()-ysig.min())/5.
        l10 = np.log10(ytick)
        if l10 < 0:
            yti = int(10**(1 + l10 - int(l10)))
#            r = -0.5
        else:
            yti = int(10**(l10 - int(l10)))
#            r = 0.5
        if yti == 3:
            yti = 2
        elif yti > 5:
            yti = 5
        ytick = yti * 10**int(np.log10(ytick/yti)+.5)
        fp.write('@{}axis tick major {}\n'.format('x',xticks[1]-xticks[0]))
        fp.write('@{}axis tick major {}\n'.format('y',ytick))
        fp.write("@type xy\n")
        for x,y,m in zip(savedX,ysig,savedX.mask):
            if not m: fp.write("{} {}\n".format(x,y))
        fp.write("&\n")
        fp.write(linedef3.format("s1",'',1,0,1.0,0,0,1))
        fp.close()
        print('file',filename,'written')

    def CopyRietveld2Origin(Pattern,Plot,Page,plotOpt,G2frame):
        # Exports plot to Origin. This function was written by Conrad Gillard (conrad.gillard@gmail.com).

        def origin_shutdown_exception_hook(exctype, value, tracebk):
            '''Ensures Origin gets shut down if an uncaught exception'''
            try: 
                op.exit()
            except:
                pass
            print('\n****OriginPro error****')
            import traceback
            traceback.print_exception(exctype, value, tracebk)
            G2G.G2MessageBox(G2frame,
                'Failed to connect to OriginPro. Is it installed?\nSee console window for more info.')
            #sys.__excepthook__(exctype, value, tracebk)

        # Function to increase line width, for later use
        def increase_line_width(plot):
            layr = op.GLayer(plot.layer)
            pindex = plot.index()
            pname = layr.obj.GetStrProp('plot{}.name'.format(pindex+1))
            layr.lt_exec('set {} -w 1000'.format(pname))
        
        #import itertools # delay this since not commonly called or needed

        try:
            import originpro as op
        except:
            note1,note2 = '',''
            # get pip location
            pyPath = os.path.split(os.path.realpath(sys.executable))[0]
            import shutil
            if shutil.which('pip'):
                pip = shutil.which('pip')
            elif os.path.exists(os.path.join(pyPath,'pip.exe')):
                pip = os.path.join(pyPath,'pip.exe')
            elif os.path.exists(os.path.join(pyPath,'Scripts','pip.exe')):
                pip = os.path.join(pyPath,'Scripts','pip.exe')
            else:
                note1 = "\nNote pip not found, you may need to install that too\n"
                pip = 'pip'
            try:
                import win32clipboard
                win32clipboard.OpenClipboard()
                win32clipboard.EmptyClipboard()
                win32clipboard.SetClipboardText('{} install originpro'.format(pip))
                win32clipboard.CloseClipboard()
                note2 = "\nNote: command copied to clipboard (use control-V to paste in cmd.exe window)"
            except:
                pass
            msg = """Use of the OriginPro exporter requires that OriginPro be
installed on your computer as well as a communication 
module (originpro) via pip. 
{}
Use command

\t{} install originpro

in a cmd.exe window to do this.
{}""".format(note1,pip,note2)
            G2G.G2MessageBox(G2frame,msg)
            return
        
        lblList = []
        valueList = []

        lblList.append('Axis-limits')
        valueList.append(list(Plot.get_xlim())+list(Plot.get_ylim()))

        tickpos = {}

        for i,l in enumerate(Plot.lines):
            lbl = l.get_label()
            if lbl[1:] in ('obs','calc','bkg','zero','diff'):
                lbl = lbl[1:]
            if 'magline' in lbl:
                pass
            elif lbl in ('obs','calc','bkg','zero','diff'):
                if lbl == 'obs':
                    lblList.append('x')
                    valueList.append(l.get_xdata())
                c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
                if sum(c) == 4.0: continue
                lblList.append(lbl)
                valueList.append(l.get_ydata())
            elif l in Page.tickDict.values():
                c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
                if sum(c) == 4.0: continue
                tickpos[lbl] = l.get_ydata()[0]
                lblList.append(lbl)
                valueList.append(l.get_xdata())
        if tickpos:
            lblList.append('tick-pos')
            valueList.append([])
            for i in tickpos:
                valueList[-1].append(i)
                valueList[-1].append(tickpos[i])
        # add (obs-calc)/sigma [=(obs-calc)*sqrt(weight)]
        lblList.append('diff/sigma')
        valueList.append(Pattern[1][5]*np.sqrt(Pattern[1][2]))
        if sum(Pattern[1][0].mask): # are there are excluded points? If so, add them too
            lblList.append('excluded')
            valueList.append(1*Pattern[1][0].mask)
        # magnifcation values
        for l in Plot.texts:
            lbl = l.get_label()
            if 'magline' not in lbl: continue
            if l.xycoords == 'axes fraction':
                lbl = 'initial-mag'
                lblList.append(lbl)
                valueList.append([l.get_text()])
            else:
                lbl = 'mag'
                lblList.append(lbl)
                valueList.append([l.get_text(),l.get_position()[0]])

        # Start Origin instance
        if op and op.oext:
            sys.excepthook = origin_shutdown_exception_hook

        # Set Origin instance visibility
        if op.oext:
            op.set_show(True)
        op.new()

        # Create folder to hold data and graph
        refinementName = G2frame.Label
        refinementName = refinementName[17:-4]
        fldr = op.po.RootFolder.Folders.Add(refinementName)
        fldr.Activate()

        # Create worksheet to hold refinement data
        dest_wks = op.new_sheet('w', lname='Refinement Data')
        dest_wks.cols = 5

        # Import refinement data
        colNamesList =["x", "Observed", 'Calculated','Background','(Obs-Calc)']
        for i in range(1, 6):
            dest_wks.from_list(col=i-1, data=valueList[i].tolist(), lname=colNamesList[i-1])

        # Create graph object, to which data will be added
        template = os.path.join(GSASIIpath.path2GSAS2,'inputs','OriginTemplate2.otpu')
        if not os.path.exists(template):  # patch 3/2024 for svn dir organization
            template = os.path.join(GSASIIpath.path2GSAS2,'OriginTemplate2.otpu')
        if not os.path.exists(template):
            print('Error: OriginTemplate2.otpu not found')
            return
        graph = op.new_graph(template=template)
        graph.lname = refinementName + "_G"
        gl = graph[0]

        # Plot observed as scatter
        plot = gl.add_plot(dest_wks, coly=1, colx=0, type='scatter')
        plot.symbol_size = 10
        plot.symbol_kind = 7
        plot.colormap = 'Classic'
        plot.color = 1

        # Plot calculated, background and difference as line
        for i in range(2, 5):
            plot = gl.add_plot(dest_wks, coly=i, colx=0, type='line')
            plot.colormap = 'Classic'
            plot.color = i
            increase_line_width(plot)

        # Import reflection data for each phase
        tickPosIdx = lblList.index("tick-pos")
        # Initialise counters for colour index and number of phases
        j = 1
        k = 1
        refLegendText = ""
        for i in range(7, tickPosIdx):
            # Create worksheet to hold reflections data
            dest_wks = op.new_sheet('w', lname=lblList[i] + " Reflections")
            dest_wks.cols = 2
            # Generate lists of tick positions
            tickPosList = valueList[i].tolist()
            # Generate lists of tick intensities
            refIntens = valueList[tickPosIdx][(i - 6) * 2 - 1]
            refIntens = float(refIntens)
            refIntensList = [refIntens] * len(tickPosList)
            # Import tick positions and intensities to worksheet
            dest_wks.from_list(col=0, data=tickPosList, lname="Tick Position")
            dest_wks.from_list(col=1, data=refIntensList, lname="Tick Intensity")
            # Add reflections to plot
            plot = gl.add_plot(dest_wks, coly=1, colx=0, type='scatter')
            plot.symbol_size = 10
            plot.symbol_kind = 10
            plot.color = 4 + j
            refLegendText = refLegendText + "\\l(" + str(4 + k) + ") " + lblList[i] + " "
            # Increment phase counter
            k += 1
            # increment colour index, skipping yellow because it cannot be seen
            if j == 2:
                j += 2
            else:
                j += 1

        #   # Set axis limits
        xmin = Plot.get_xlim()[0]
        if Plot.dataLim.x0 > xmin:
            xmin = Plot.dataLim.x0 + 0.1
        xmax = Plot.get_xlim()[1]
        if Plot.dataLim.x1 < xmax:
            xmax = Plot.dataLim.x1 + 0.1
        gl.set_xlim(xmin, xmax)

        ymin = Plot.get_ylim()[0]
        ymax = Plot.get_ylim()[1]
        gl.set_ylim(ymin, ymax)

        # Change graph titles
        gl.axis(ax="x").title = "2θ (°)"
        gl.axis(ax="y").title = "Intensity (Arbitrary units)"

        # Set up legend
        label = gl.label('Legend')
        label.text = '\\l(1) %(1)\\l(2) %(2)\\l(3) %(3)\\l(4) %(4) %(CRLF)' + refLegendText

    def CopyRietveld2Igor(Pattern,Plot,Page,plotOpt,filename,G2frame):
        '''Copy the contents of the Rietveld graph from the plot window to
        a Igor Pro input file (tested with v7.05). Uses values from Pattern 
        to also generate a delta/sigma plot below. 

        Coded with lots of help from Jan Ilavsky.
        '''
        import itertools # delay this until called, since not commonly needed
        InameDict = {'obs':'Intensity', 'calc':'FitIntensity',
                        'bkg':'Background','diff':'Difference',
                        'omcos':'NormResidual','zero':'Zero'}
        igor_symbols = {"o":19, "s":18, "D":29, "^":17, "3":46, 'v':23,
                "4":49, "+":1, "P":60, "x":0, "X":62, "*":2}
        def Write2cols(fil,dataItems):
            '''Write a line to a file in space-separated columns. 
            Skips masked items. 

            :param object fil: file object
            :param list dataItems: items to write as row in file
            '''
            line = ''
            for item in dataItems:
                if ma.is_masked(item): return
                if line: line += ' '
                item = str(item)
                if ' ' in item:
                    line += '"'+item+'"'
                else:
                    line += item
            fil.write(line+'\n')
        proj = os.path.splitext(G2frame.GSASprojectfile)[0]
        if not proj: proj = 'GSASIIproject'
        proj = proj.replace(' ','')
        valueList = []
        markerSettings = []
        Icolor = {}
        legends = []
        zerovals = None
        fontsize = 18*float(plotOpt['labelSize'])/12.
        for i,l in enumerate(Plot.lines):
            lbl = l.get_label()
            if not plotOpt['Show'].get(lbl[1:],True): continue
            if plotOpt['legend'].get(lbl[1:]):
                legends.append((InameDict[lbl[1:]],lbl[1:]))
            if l.get_label()[1:] in ('obs','calc','bkg','zero','diff'):
                lbl = l.get_label()[1:]
            if plotOpt['legend'].get(lbl) and lbl in InameDict:
                legends.append((InameDict[lbl],lbl))
            if lbl == 'obs':
                x = l.get_xdata()
                valueList.append(x)
                zerovals = (x.min(),x.max())
                gsiz = 5*float(plotOpt['markerSiz'])/8.
                marker = plotOpt['markerSym']
                gsym = igor_symbols.get(marker,12)
                gmw = float(plotOpt['markerWid'])
                markerSettings.append(
                    'mode({0})=3,marker({0})={1},msize({0})={2},mrkThick({0})={3}'
                    .format('Intensity',gsym,gsiz,gmw))
            elif lbl in ('calc','bkg','zero','diff'):
                markerSettings.append(
                    'mode({0})=0, lsize({0})={1}'
                    .format(InameDict[lbl],plotOpt['lineWid']))
            else:
                continue
            c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
            Icolor[InameDict[lbl]] = [j*65535 for j in c[0:3]]
            if lbl != 'zero':
                valueList.append(l.get_ydata())
        valueList.append(Pattern[1][5]*np.sqrt(Pattern[1][2]))
        fp = open(filename,'w')
        fp.write('''IGOR
X setDataFolder root:
X //   ***   Replace GSAS2Data with name of current project in GSAS. 
X //   ***   this name will get "order" number (0,1,2...) to be unique and data will be stored there. 
X //   ***   and the graph will also be named using this base name. 
X NewDataFolder/O/S $(UniqueName(CleanupName("{}",0),11, 0))
X string GSAXSProjectName = GetDataFolder(0)
WAVES /D/O TwoTheta, Intensity, FitIntensity, Background, Difference, NormResidual
BEGIN
'''.format(proj))
        # invert lists into columns, use iterator for all values
        for row in itertools.zip_longest(*valueList,fillvalue=' '):
            Write2cols(fp,row)
        fp.write('''END
X //  ***   static part of the code, NB reflection tickmarks later ****
X SetScale d 0,0, "degree", TwoTheta
X //  ***   this is where graph is created and data added
X string G_Name=CleanupName(GSAXSProjectName,0)
X Display/K=1/W=(50,40,850,640) Intensity vs twoTheta; DoWindow/C $(G_Name) 
X DoWindow/T $(G_Name), G_Name
X AppendToGraph FitIntensity vs TwoTheta
X AppendToGraph Background vs TwoTheta
X AppendToGraph Difference vs TwoTheta
X AppendToGraph/L=Res_left/B=Res_bot NormResidual vs TwoTheta
X //  ***   Here we have modification of the four axes used in the graph
X ModifyGraph mode=2,lsize=5, mirror=2
X SetAxis/A/E=2 Res_left
X ModifyGraph freePos(Res_left)={0,kwFraction}
X ModifyGraph freePos(Res_bot)={0.2,kwFraction}
X ModifyGraph standoff(bottom)=0
X ModifyGraph axisEnab(left)={0.2,1}
X ModifyGraph axisEnab(Res_left)={0,0.2}
X ModifyGraph lblPosMode(Res_left)=1, mirror(Res_bot)=0
X ModifyGraph tickUnit(bottom)=1,manTick=0,manMinor(bottom)={0,50}
X ModifyGraph tick(Res_bot)=1,noLabel(Res_bot)=2,standoff(Res_bot)=0
X ModifyGraph manMinor(Res_bot)={0,50},tick(Res_bot)=2
X ModifyGraph axThick=2
X ModifyGraph mirror(Res_bot)=0,nticks(Res_left)=3,highTrip(left)=1e+06,manTick(bottom)=0
X ModifyGraph btLen=5
''')
        fp.write('X ModifyGraph gfSize={}\n'.format(int(fontsize+.5)))
        
        # line at zero
        if not zerovals:
            zerovals = (Plot.get_xlim()[0],Plot.get_xlim()[1])
        fp.write('X //  ***   add line at y=zero\n')
        fp.write('WAVES /D/O ZeroX, Zero\nBEGIN\n')
        fp.write(' {0} 0.0\n {1} 0.0\n'.format(*zerovals))
        fp.write('END\nX AppendToGraph Zero vs ZeroX\n')
        if 'sqrt' in Plot.yaxis.get_label().get_text():
            ylabel = '\u221AIntensity'            
        else:
            ylabel = 'Intensity'
        fp.write('''X //  ***   add axis labels and position them
X Label left "{1}"
X Label Res_left "{2}"
X Label bottom "{0}"
X ModifyGraph lblPosMode=0,lblPos(Res_left)=84
'''.format("2Θ",ylabel,"∆/σ"))
        fp.write('''X //  ***   set display limits. 
X SetAxis left {2}, {3}
X SetAxis bottom {0}, {1}
X SetAxis Res_bot {0}, {1}
'''
                    .format(Plot.get_xlim()[0],Plot.get_xlim()[1],
                                    Plot.get_ylim()[0],Plot.get_ylim()[1]))
        fp.write('X //  ***   modify how data are displayed ****\n')
#         fp.write(
# '''X //  ***   here starts mostly static part of the code, here we modify how data are displayed ****
# X ModifyGraph mode(FitIntensity)=0,rgb(FitIntensity)=(1,39321,19939), lsize(FitIntensity)=1
# X ModifyGraph mode(Background)=0,lsize(Background)=2, rgb(Background)=(65535,0,0)
# X ModifyGraph mode(Difference)=0,lsize(Difference)=2,rgb(Difference)=(0,65535,65535)
# X //  ***   modifications for the bottom graph, here we modify how data are displayed ****
# X ModifyGraph mode(NormResidual)=0,lsize(NormResidual)=1,rgb(NormResidual)=(0,0,0)
# X //  ***   end of modifications for the main data in graph
# ''')
        for m in markerSettings:
                fp.write('X ModifyGraph {}\n'.format(m))
        for lbl in Icolor:
                fp.write('X ModifyGraph rgb({})=({:.0f},{:.0f},{:.0f})\n'
                                 .format(lbl,*Icolor[lbl]))
        fp.write('X ModifyGraph mode(NormResidual)=0,lsize(NormResidual)=1,rgb(NormResidual)=(0,0,0)\n')
        fp.write('X //  ***   End modify how data are displayed ****\n')
        # loop over reflections
        ticknum = 0
        for i,l in enumerate(Plot.lines):
            lbl = l.get_label()
            if not plotOpt['Show'].get(lbl,True): continue
            if l in Page.tickDict.values():
                c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
                if sum(c) == 4.0: continue # white is invisible
                ticknum += 1
                phasename = 'tick{}'.format(ticknum)
                if plotOpt['legend'].get(lbl):
                    legends.append((phasename,plotOpt['phaseLabels'].get(lbl,lbl)))
                tickpos = l.get_ydata()[0]
                fp.write(
'''X //  reflection tickmark for phase {1}
WAVES /D/O {0}X, {0}
BEGIN
'''.format(phasename,lbl))
                for x in l.get_xdata():
                    fp.write('{} {}\n'.format(x,tickpos))
                fp.write('''END
X //  ***   controls for {1} reflection tickmarks
X AppendToGraph {0} vs {0}X
X ModifyGraph mode({0})=3,mrkThick({0})=2,gaps({0})=0
X ModifyGraph marker({0})=10,rgb({0})=({2},{3},{4})
'''.format(phasename,lbl,*[j*65535 for j in c[0:3]]))
                
        # plot magnification lines and labels
        j = 0
        for i,l in enumerate(Plot.lines):
            lbl = l.get_label()
            if 'magline' not in lbl: continue
            j += 1
            fp.write('WAVES /D/O mag{0}X, mag{0}\nBEGIN\n'.format(j))
            fp.write(' {0} {1}\n {0} {2}\n'.format(
                l.get_data()[0][0],Plot.get_ylim()[0],Plot.get_ylim()[1]))
            fp.write('END\nX AppendToGraph mag{0} vs mag{0}X\n'.format(j))
            fp.write('X ModifyGraph lstyle(mag{0})=3,rgb(mag{0})=(0,0,0)\n'.format(j))
        for l in Plot.texts:
            if 'magline' not in l.get_label(): continue
            if l.xycoords[0] == 'data':
                xpos = l.get_position()[0]
            elif l.get_position()[0] == 0:
                xpos = Plot.get_xlim()[0]
            else:
                xpos = Plot.get_xlim()[1]*.95
            fp.write('X SetDrawEnv xcoord=bottom, ycoord= abs,textyjust=2,fsize={}\n'.format(fontsize))
            fp.write('X DrawText {0},2,"{1}"\n'.format(xpos,l.get_text()))

        # legend
        s = ""
        for nam,txt in legends:
            if s: s += r'\r'
            s += r'\\s({}) {}'.format(nam,txt)
        fp.write('X Legend/C/N=text0/J "{}"\n'.format(s))
        fp.close()
        print('file',filename,'written')

    def CopyRietveld2csv(Pattern,Plot,Page,filename):
        '''Copy the contents of the Rietveld graph from the plot window to
        .csv file
        '''
        import itertools # delay this since not commonly called or needed

        lblList = []
        valueList = []
        tickpos = {}
        lblList.append('used')
        try:
            valueList.append([0 if i else 1 for i in savedX.mask])
        except:
            pass
        if 'TOF' in Plot.xaxis.get_label_text():
            lblList.append('x, TOF (msec)')
        elif 'Q,' in Plot.xaxis.get_label_text():
            lblList.append('x, Q (A-1)')
        elif 'd,' in Plot.xaxis.get_label_text():
            lblList.append('x, d-space (A)')
        elif 'E,' in Plot.xaxis.get_label_text():
            lblList.append('x, E (keV)')
        elif 'theta' in Plot.xaxis.get_label_text():
            lblList.append('x, 2theta (deg)')
        else:
            lblList.append('x, ?')
        valueList.append(savedX.data)
        for i,l in enumerate(Plot.lines):
            lbl = l.get_label()
            if lbl[1:] in ('obs','calc','bkg','zero','diff'):
                lbl = l.get_label()[1:]
            if 'magline' in lbl:
                pass
            elif lbl in ('obs','calc','bkg','zero','diff'):
                c = plotOpt['colors'].get(lbl,l.get_color())
                if sum(c) == 4.0: continue  #skip over any "white" entries
                if lbl == 'zero': continue
                if lbl == 'obs':  # include all observed data
                    lblList.append('obs')
                    valueList.append(Pattern[1][1].data)
                    continue
                lblList.append(lbl)
                valueList.append(l.get_ydata())
            elif l in Page.tickDict.values():
                c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
                if sum(c) == 4.0: continue
                tickpos[lbl] = l.get_ydata()[0]
                lblList.append(lbl)
                valueList.append(l.get_xdata())
        if tickpos:
            lblList.append('tick-pos')
            valueList.append([])
            for i in tickpos:
                valueList[-1].append(i)
                valueList[-1].append(tickpos[i])
        # add (obs-calc)/sigma [=(obs-calc)*sqrt(weight)]
        lblList.append('diff/sigma')
        wtFactor = Pattern[0]['wtFactor']
        DZ = (Pattern[1][1]-Pattern[1][3])*np.sqrt(wtFactor*Pattern[1][2])
        valueList.append(DZ)

        if sum(Pattern[1][0].mask): # are there are excluded points? If so, add them too
            lblList.append('excluded')
            valueList.append(1*Pattern[1][0].mask)
        lblList.append('Axis-limits')
        valueList.append(list(Plot.get_xlim())+list(Plot.get_ylim()))
        # magnification values
        for l in Plot.texts:
            lbl = l.get_label()
            if 'magline' not in lbl: continue
            if l.xycoords == 'axes fraction':
                lbl = 'initial-mag'
                lblList.append(lbl)
                valueList.append([l.get_text()])
            else:
                lbl = 'mag'
                lblList.append(lbl)
                valueList.append([l.get_text(),l.get_position()[0]])
        fp = open(filename,'w')
        G2plt.Write2csv(fp,lblList,header=True)
        # invert lists into columns, use iterator for all values
        for row in itertools.zip_longest(*valueList,fillvalue=' '):
            G2plt.Write2csv(fp,row)
        fp.close()
        
    def onSave(event):
        '''Write the current plot to a file
        '''
        hcfigure = mplfig.Figure(dpi=plotOpt['dpi'],figsize=(plotOpt['width'],plotOpt['height']))
        CopyRietveldPlot(G2frame,Pattern,Plot,Page,hcfigure,plotWidgets['phaseList'])
        if 'OriginPro' in plotOpt['format']:
            CopyRietveld2Origin(Pattern,Plot,Page,plotOpt,G2frame)
            dlg.EndModal(wx.ID_OK)
            return
        longFormatName,typ = plotOpt['format'].split(',')
        fil = G2G.askSaveFile(G2frame,'','.'+typ.strip(),longFormatName)
        if 'csv' in typ and fil:
            CopyRietveld2csv(Pattern,Plot,Page,fil)
        elif 'agr' in typ and fil:
            CopyRietveld2Grace(Pattern,Plot,Page,plotOpt,fil)
        elif 'itx' in typ and fil:
            CopyRietveld2Igor(Pattern,Plot,Page,plotOpt,fil,G2frame)
        elif fil:
            if hcfigure.canvas is None:
                if GSASIIpath.GetConfigValue('debug'): print('creating canvas')
                hcCanvas(hcfigure)
            hcfigure.savefig(fil,format=typ.strip())
        dlg.EndModal(wx.ID_OK)
            
    def OnSelectColour(event):
        '''Respond to a change in color
        '''
#        lbl = plotOpt['colorButtons'].get(list(event.GetEventObject())[:3])
        if event.GetEventObject() not in plotWidgets['colorButtons']:
            print('Unexpected button',str(event.GetEventObject()))
            return
        lbl = plotWidgets['colorButtons'][event.GetEventObject()]
        c = event.GetValue()
        plotOpt['colors'][lbl] = (c.Red()/255.,c.Green()/255.,c.Blue()/255.,c.alpha/255.)
        RefreshPlot()

    def OnSavePlotOpt(event):
        'Save current "publish" settings to a JSON file with extension .pltopts'
        import json
        filroot = os.path.splitext(G2frame.GSASprojectfile)[0]+'.pltopts'
        saveFile = G2G.askSaveFile(G2frame,filroot,'.pltopts',
                                        'Saved plot options')
        if saveFile:
            json.dump(plotOpt,open(saveFile,'w'))
            print(f'Plot options written to {saveFile!r}')
            
    def OnLoadPlotOpt(event):
        'Set current "publish" settings from saved settings in a JSON file (extension .pltopts)'
        import json
        filroot = os.path.splitext(G2frame.GSASprojectfile)[0]+'.pltopts'
        fdlg = wx.FileDialog(dlg, 'Choose saved plot options file',
                    filroot,
                    style=wx.FD_OPEN, wildcard='plot options file(*.pltopts)|*.pltopts')
        if fdlg.ShowModal() == wx.ID_OK:
            filroot = fdlg.GetPath()
            plotOpt.update(json.load(open(filroot,'r')))
            wx.CallAfter(PublishRietveldPlot,G2frame,Pattern,Plot,Page,reuse=dlg)
        fdlg.Destroy()

    #### start of PublishRietveldPlot
    plotWidgets = {}
    plotWidgets['phaseList'] = []
    Init_fmts()
    if reuse:
        dlg = reuse
        vbox = dlg.GetSizer()
        vbox.Clear()
    else:
        if plotOpt['initNeeded']: Initialize()
        GetColors()            
        dlg = wx.Dialog(G2frame.plotFrame,title="Publication plot creation",
                style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        vbox = wx.BoxSizer(wx.VERTICAL)
    
    # size choices
    symChoices = ('+','x','.','o','^','v','*','|')
    txtChoices = [str(i) for i in range (8,26)]
    sizChoices = [str(i) for i in range (2,21)]
    lwidChoices = ('0.5','0.7','1','1.5','2','2.5','3','4')
    sizebox = wx.BoxSizer(wx.HORIZONTAL)
    sizebox.Add(wx.StaticText(dlg,wx.ID_ANY,'Text size'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,txtChoices,None,None,plotOpt,'labelSize',RefreshPlot,
                                   size=(50,-1))
    sizebox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    sizebox.Add((1,1),1,wx.EXPAND,1)
    sizebox.Add(wx.StaticText(dlg,wx.ID_ANY,' Obs type'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,symChoices,None,None,plotOpt,'markerSym',RefreshPlot,
                                   size=(40,-1))
    sizebox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    sizebox.Add(wx.StaticText(dlg,wx.ID_ANY,' size'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,sizChoices,None,None,plotOpt,'markerSiz',RefreshPlot,
                                   size=(50,-1))
    sizebox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    sizebox.Add(wx.StaticText(dlg,wx.ID_ANY,' width'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,lwidChoices,None,None,plotOpt,'markerWid',RefreshPlot,
            size=(50,-1))
    sizebox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    sizebox.Add((1,1),1,wx.EXPAND,1)
    sizebox.Add(wx.StaticText(dlg,wx.ID_ANY,' Line widths'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,lwidChoices,None,None,plotOpt,'lineWid',RefreshPlot,
            size=(50,-1))
    sizebox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    sizebox.Add((1,1),1,wx.EXPAND,1)
    sizebox.Add(wx.StaticText(dlg,wx.ID_ANY,' Tick size'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,sizChoices,None,None,plotOpt,'tickSiz',RefreshPlot,
            size=(50,-1))
    sizebox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    sizebox.Add(wx.StaticText(dlg,wx.ID_ANY,' width'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,lwidChoices,None,None,plotOpt,'tickWid',RefreshPlot,
            size=(50,-1))
    sizebox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    sizebox.Add((1,1),1,wx.EXPAND,1)
    helpinfo = '''----   Help on creating hard copy   ----
    Select options such as the size of text and colors for plot contents here.
    
    Tricks:
    * Use a color of pure white to remove an element from the plot (light
    gray is plotted)
    * LaTeX-like coding can be used for phase labels such as
    $\\rm FeO_2$ (for a subscript 2) or $\\gamma$-Ti for a Greek "gamma"
    
    Note that the dpi value is ignored for svg and pdf files, which are
    drawn with vector graphics (infinite resolution). Likewise, the agr and 
    itx options create input files for programs Grace (QtGrace) and Igor 
    Pro, that largely duplicate the displayed plot. Once read into Grace
    or Igor the graphs can then be customized.
    '''
    if sys.platform == "win32":
        helpinfo += '''
Note that the OriginPro connection export requires Origin 2021 or later.'''
    hlp = G2G.HelpButton(dlg,helpinfo)
    sizebox.Add(hlp,0,wx.ALL)
    vbox.Add(sizebox,0,wx.ALL|wx.EXPAND)
    
    # table of colors and legend options
    cols = 1+len(plotOpt['lineList']) + len(plotWidgets['phaseList'] )
    import wx.lib.scrolledpanel as wxscroll
    gpanel = wxscroll.ScrolledPanel(dlg,size=(600,105))
    gsizer = wx.FlexGridSizer(cols=cols,hgap=2,vgap=2)
    gsizer.Add((-1,-1))
    for lbl in plotOpt['lineList']:
        gsizer.Add(wx.StaticText(gpanel,wx.ID_ANY,lbl),0,wx.ALL)
    for lbl in plotWidgets['phaseList']:
        if lbl not in plotOpt['phaseLabels']: plotOpt['phaseLabels'][lbl] = lbl
        val = G2G.ValidatedTxtCtrl(gpanel,plotOpt['phaseLabels'],lbl,size=(110,-1),
                                   style=wx.TE_CENTRE,OnLeave=RefreshPlot)
        gsizer.Add(val,0,wx.ALL)
    gsizer.Add(wx.StaticText(gpanel,wx.ID_ANY,'Show'),0,wx.ALL)
    for lbl in list(plotOpt['lineList']) + list(plotWidgets['phaseList'] ):
        if lbl not in plotOpt['Show']:
            plotOpt['Show'][lbl] = True
        ch = G2G.G2CheckBox(gpanel,'',plotOpt['Show'],lbl,RefreshPlot)
        gsizer.Add(ch,0,wx.ALL|wx.ALIGN_CENTER)
    gsizer.Add(wx.StaticText(gpanel,wx.ID_ANY,'Include in legend'),0,wx.ALL)
    for lbl in list(plotOpt['lineList']) + list(plotWidgets['phaseList'] ):
        if lbl not in plotOpt['legend']:
            plotOpt['legend'][lbl] = False
        ch = G2G.G2CheckBox(gpanel,'',plotOpt['legend'],lbl,RefreshPlot)
        gsizer.Add(ch,0,wx.ALL|wx.ALIGN_CENTER)
    gsizer.Add(wx.StaticText(gpanel,wx.ID_ANY,'Color'),0,wx.ALL)
    plotWidgets['colorButtons'] = {}
    for lbl in list(plotOpt['lineList']) + list(plotWidgets['phaseList']):
        import  wx.lib.colourselect as csel
        if lbl not in  plotOpt['colors']:
            plotOpt['colors'][lbl] = (0.5, 0.5, 0.5, 1)
        color = wx.Colour(*[int(255*i) for i in plotOpt['colors'][lbl]])
        b = csel.ColourSelect(gpanel, -1, '', color)
        b.Bind(csel.EVT_COLOURSELECT, OnSelectColour)
        plotWidgets['colorButtons'][b] = lbl
        gsizer.Add(b,0,wx.ALL|wx.ALIGN_CENTER)
    # plot limits
    vlbox = wx.BoxSizer(wx.VERTICAL)
    vlbox.Add(wx.StaticText(dlg,wx.ID_ANY,'plot limit overrides'),0,wx.TOP|wx.ALIGN_CENTER,5)
    limsizer = wx.FlexGridSizer(cols=3,hgap=2,vgap=2)
    for xy in 'x','y':
        for minmax in 'min','max':
            key = f'{xy}{minmax}'
            txt = wx.StaticText(dlg,wx.ID_ANY,key)
            limsizer.Add(txt,0,wx.ALIGN_RIGHT)
            ch = G2G.G2CheckBox(dlg,'',plotOpt,key+'_use',RefreshPlot)
            limsizer.Add(ch,0,wx.ALL|wx.ALIGN_CENTER)
            val = G2G.ValidatedTxtCtrl(dlg,plotOpt,key,size=(90,-1),
                                   style=wx.TE_CENTRE,OnLeave=RefreshPlot)
            limsizer.Add(val,0,wx.ALL)
    vlbox.Add(limsizer,0,wx.BOTTOM,5)
    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add((5,-1))
    hbox.Add(vlbox,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    # assemble controls box
    hbox.Add((10,-1))
    hbox.Add((-1,-1),1,wx.EXPAND,0)
    gpanel.SetSizer(gsizer)
    gpanel.SetAutoLayout(1)
    gpanel.SetupScrolling()
    hbox.Add(gpanel,0,wx.ALIGN_CENTER)
    hbox.Add((-1,-1),1,wx.EXPAND,0)
    vbox.Add(hbox,0,wx.ALL|wx.EXPAND)

    # hard copy options
    hbox = wx.BoxSizer(wx.HORIZONTAL)
    txt = wx.StaticText(dlg,wx.ID_ANY,'Hard copy  ')
    hbox.Add(txt,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    txt = wx.StaticText(dlg,wx.ID_ANY,' Width (in):')
    hbox.Add(txt,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    val = G2G.ValidatedTxtCtrl(dlg,plotOpt,'width',xmin=3.,xmax=20.,nDig=(5,1),size=(45,-1))
    hbox.Add(val,0,wx.ALL)
    txt = wx.StaticText(dlg,wx.ID_ANY,' Height (in):')
    hbox.Add(txt,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    val = G2G.ValidatedTxtCtrl(dlg,plotOpt,'height',xmin=3.,xmax=20.,nDig=(5,1),size=(45,-1))
    hbox.Add(val,0,wx.ALL)
    txt = wx.StaticText(dlg,wx.ID_ANY,'File format:')
    hbox.Add(txt,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    fmtval = G2G.EnumSelector(dlg,plotOpt,'format',plotWidgets['fmtChoices'])
    hbox.Add(fmtval,0,wx.ALL)
    # TODO: would be nice to gray this out when format will ignore this
    ptxt = wx.StaticText(dlg,wx.ID_ANY,'Pixels/inch:')
    pval = G2G.ValidatedTxtCtrl(dlg,plotOpt,'dpi',xmin=60,xmax=1600,
                                    size=(40,-1))
    hbox.Add(ptxt,0,wx.ALL|wx.ALIGN_CENTER_VERTICAL)
    hbox.Add(pval,0,wx.ALL)
    vbox.Add(hbox,0,wx.ALL|wx.ALIGN_CENTER)

    # screen preview
    figure = mplfig.Figure(figsize=(plotOpt['width'],plotOpt['height']))
    canvas = Canvas(dlg,-1,figure)
    vbox.Add(canvas,1,wx.ALL|wx.EXPAND,1)

    # buttons at bottom
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    btnsizer.Add((5,-1))
    btn = wx.Button(dlg, wx.ID_ANY,'save\nsettings',size=(80,-1))
    btn.Bind(wx.EVT_BUTTON,OnSavePlotOpt)
    btnsizer.Add(btn)
    btnsizer.Add((5,-1))
    btn = wx.Button(dlg, wx.ID_ANY,'load\nsettings',size=(80,-1))
    btn.Bind(wx.EVT_BUTTON,OnLoadPlotOpt)
    btnsizer.Add(btn)
    btnsizer.Add((-1,-1),1,wx.EXPAND)
    btn = wx.Button(dlg, wx.ID_CANCEL)
    btnsizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
    btnsizer.Add((5,-1))
    btn = wx.Button(dlg, wx.ID_SAVE)
    btn.Bind(wx.EVT_BUTTON,onSave)
    btnsizer.Add(btn,0,wx.ALIGN_CENTER_VERTICAL)
    btnsizer.Add((-1,-1),1,wx.EXPAND)
    btnsizer.Add((170,-1),1,wx.EXPAND) # empty box to center buttons
    vbox.Add(btnsizer, 0, wx.EXPAND|wx.TOP|wx.BOTTOM,4)
    dlg.SetSizer(vbox)
    vbox.Fit(dlg)
    dlg.Layout()
    dlg.CenterOnParent()

    CopyRietveldPlot(G2frame,Pattern,Plot,Page,figure,plotWidgets['phaseList'])     # preview plot
    figure.canvas.draw()

    if not reuse:
        dlg.ShowModal()
        dlg.Destroy()
    fmtval.SetFocus()   # move focus off from pixel size
    return

def CopyRietveldPlot(G2frame,Pattern,Plot,Page,figure,phaseList):
    '''Copy the contents of the Rietveld graph from the plot window to another
    mpl figure which can be on screen or can be a file for hard copy.
    Uses values from Pattern to also generate a delta/sigma plot below the 
    main figure, since the weights are not available from the plot. 

    :param list Pattern: histogram object from data tree
    :param mpl.axes Plot: The axes object from the Rietveld plot
    :param wx.Panel Page: The tabbed panel for the Rietveld plot
    :param matplotlib.figure.Figure figure: The figure object from the Rietveld plot
    '''
    # set up axes
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax0 = figure.add_subplot(gs[0])
    ax1 = figure.add_subplot(gs[1])
    figure.subplots_adjust(left=int(plotOpt['labelSize'])/100.,bottom=int(plotOpt['labelSize'])/150.,
                           right=.98,top=1.-int(plotOpt['labelSize'])/200.,hspace=0.0)
    ax0.tick_params('x',direction='in',labelbottom=False)
    ax0.tick_params(labelsize=plotOpt['labelSize'])
    ax1.tick_params(labelsize=plotOpt['labelSize'])
    if mpl.__version__.split('.')[0] == '1': # deal with older matplotlib, which puts too many ticks
        ax1.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=2))
        ax1.yaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins=4))
    ax1.set_xlabel(Plot.get_xlabel(),fontsize=plotOpt['labelSize'])
    ax0.set_ylabel(Plot.get_ylabel(),fontsize=plotOpt['labelSize'])
    ax1.set_ylabel(r'$\Delta/\sigma$',fontsize=plotOpt['labelSize'])
    # set axes ranges, get ranges from display
    lims = {}
    lims['xmin'],lims['xmax'] = Plot.get_xlim()
    lims['ymin'],lims['ymax'] = Plot.get_ylim()
    for key in 'xmin','xmax','ymin','ymax':   # apply limit overrides where use flag is set
        if plotOpt[key+'_use']: lims[key] = plotOpt[key]
    ax0.set_xlim((lims['xmin'],lims['xmax']))
    ax1.set_xlim((lims['xmin'],lims['xmax']))
    ax0.set_ylim((lims['ymin'],lims['ymax']))
    
    legLbl = []
    legLine = []
    # get the obs/calc... & magnification lines and xfer them
    for i,l in enumerate(Plot.lines):
        lbl = l.get_label()
        if lbl[1:] in ('obs','calc','bkg','zero','diff'):
            lbl = lbl[1:]
        if 'magline' in lbl:                              # magnification lines
            ax0.axvline(l.get_data()[0][0],color='0.5',dashes=(1,1))
        elif lbl in ('obs','calc','bkg','zero','diff'):   # data segments
            if not plotOpt['Show'].get(lbl,True): continue
            marker = l.get_marker()
            lineWid = l.get_lw()
            siz = l.get_markersize()
            mew = l.get_mew()
            if lbl == 'obs':
                siz = float(plotOpt['markerSiz'])
                marker = plotOpt['markerSym']
                mew = float(plotOpt['markerWid'])
            else:
                lineWid = float(plotOpt['lineWid'])
            c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
            if sum(c) == 4.0: continue
            if plotOpt['legend'].get(lbl):
                uselbl = lbl
            else:
                uselbl = '_'+lbl
            if lbl == 'zero':
                art = [ax0.axhline(0.,color=c,
                     lw=lineWid,label=uselbl,ls=l.get_ls(),
                     marker=marker,ms=siz,mew=mew)]
            else:
                art = ax0.plot(l.get_xdata(),l.get_ydata(),color=c,
                     lw=lineWid,label=uselbl,ls=l.get_ls(),
                     marker=marker,ms=siz,mew=mew,
                     )
            if plotOpt['legend'].get(lbl):
                legLbl.append(uselbl)
                legLine.append(art[0])
        elif l in Page.tickDict.values():                 # reflection tickmarks
            if not plotOpt['Show'].get(lbl,True): continue
            c = plotOpt['colors'].get(lbl,mpl.colors.to_rgba(l.get_color()))
            #siz = l.get_markersize()
            siz = float(plotOpt['tickSiz'])
            #mew = l.get_mew()
            mew = float(plotOpt['tickWid'])
            if sum(c) == 4.0: continue
            if not plotOpt['legend'].get(lbl):
                uselbl = '_'+lbl
            else:
                uselbl = plotOpt['phaseLabels'].get(lbl,lbl)
            art = ax0.plot(l.get_xdata(),l.get_ydata(),color=c,
                     lw=l.get_lw(),ls=l.get_ls(),label=uselbl,
                     marker=l.get_marker(),ms=siz,mew=mew,
                     )
            if plotOpt['legend'].get(lbl):
                legLbl.append(uselbl)
                legLine.append(art[0])
        elif '_FLT_' in lbl:                              # full-length vertical reflection markers
            ph = lbl[5:]
            c = plotOpt['colors'][ph]
            ax0.axvline(l.get_xdata()[0],color=c,lw=float(plotOpt['tickWid']))
        elif lbl in plotOpt['colors']:                    # draw phase partials
            if l.get_marker() != 'None': continue # ignore 2nd tickmark artist
            if not plotOpt['Show'].get(lbl,True): continue # remove partial w/tickmarks
            c = plotOpt['colors'][lbl]
            lineWid = float(plotOpt['lineWid'])
            # xfer the phase partials to the new plot, overriding the phase color and line width
            ax0.plot(l.get_xdata(),l.get_ydata(),color=c,
                linewidth=lineWid,
                linestyle=l.get_ls())
#        elif GSASIIpath.GetConfigValue('debug'):
#            print('other line:',lbl)

    # copy text items: magnification labels and reflection markers
    for l in Plot.texts:
        lbl = l.get_label()
        if lbl[1:] in phaseList:                           # reflection markers
            if not plotOpt['Show'].get(lbl[1:],True): continue # remove marker w/tickmarks
            p = l.get_bbox_patch()
            props = {'facecolor':p.get_facecolor(),
                     'edgecolor':p.get_facecolor(),
                     'alpha':p.get_alpha(),
                     'boxstyle':p.get_boxstyle()}
            c = plotOpt['colors'][lbl[1:]]
            ax0.text(l.get_position()[0],l.get_position()[1],
                      l.get_text(),
                      fontsize=float(plotOpt['labelSize']),
                      c=c,
                      ha='center',
                      va='top',
                      bbox=props,
                      rotation=l.get_rotation())
        elif 'maglbl' in lbl:                             # label on a magnification region
            ax0.annotate(l.get_text(),
                    xy=(l.get_position()), xycoords=l.xycoords,
                    verticalalignment='bottom',
                    horizontalalignment=l.get_horizontalalignment(),
                    fontsize=float(plotOpt['labelSize']))
#        elif GSASIIpath.GetConfigValue('debug'):
#            print('other text:',l.get_label())
            
    # generate the (obs-calc)/sigma values and plot them
    wtFactor = Pattern[0]['wtFactor']
    DZ = (Pattern[1][1]-Pattern[1][3])*np.sqrt(wtFactor*Pattern[1][2])
    ax1.plot(savedX,DZ,color='k')
    # show the legend, if anything is in it (list legLine)
    if legLine:
        ax0.legend(legLine,legLbl,loc='best',prop={'size':plotOpt['labelSize']})
    
def changePlotSettings(G2frame,Plot):
    '''Code in development to allow changes to plot settings
    prior to export of plot with "floppy disk" button
    '''
    def RefreshPlot(*args,**kwargs):
        '''Apply settings to the plot
        '''
        Plot.figure.subplots_adjust(left=int(plotOpt['labelSize'])/100.,
                            bottom=int(plotOpt['labelSize'])/150.,
                            right=.98,
                            top=1.-int(plotOpt['labelSize'])/200.,
                            hspace=0.0)
        for P in Plot.figure.axes:
            P.get_xaxis().get_label().set_fontsize(plotOpt['labelSize'])
            P.get_yaxis().get_label().set_fontsize(plotOpt['labelSize'])
            for l in P.get_xaxis().get_ticklabels():
                l.set_fontsize(plotOpt['labelSize'])
            for l in P.get_yaxis().get_ticklabels():
                l.set_fontsize(plotOpt['labelSize'])
            for l in P.lines: 
                l.set_linewidth(float(plotOpt['lineWid']))
            P.get_xaxis().set_tick_params(width=float(plotOpt['lineWid']))
            P.get_yaxis().set_tick_params(width=float(plotOpt['lineWid']))
            for l in P.spines.values():
                l.set_linewidth(float(plotOpt['lineWid']))
        if Page.plotStyle.get('title',True):        
            Plot.set_title(plotOpt['title'],fontsize=plotOpt['labelSize'])
        for i,P in enumerate(Plot.figure.axes):
            if not P.get_visible(): continue
            if i == 0:
                lbl = ''
            else:
                lbl = str(i)
            P.get_xaxis().set_label_text(plotOpt['xtitle'+lbl])
            P.get_yaxis().set_label_text(plotOpt['ytitle'+lbl])            
        Plot.figure.canvas.draw()

    txtChoices = [str(i) for i in range (8,26)]
    lwidChoices = ('0.5','0.7','1','1.5','2','2.5','3','4')
    dlg = wx.Dialog(G2frame.plotFrame,style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
    vbox = wx.BoxSizer(wx.VERTICAL)
    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add(wx.StaticText(dlg,wx.ID_ANY,'Text size'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,txtChoices,None,None,plotOpt,'labelSize',RefreshPlot,
                                   size=(50,-1))
    hbox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    vbox.Add(hbox,0,wx.ALL|wx.EXPAND)
    
    vbox.Add((1,5))
    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add(wx.StaticText(dlg,wx.ID_ANY,' Line widths'),0,wx.ALL)
    w = G2G.G2ChoiceButton(dlg,lwidChoices,None,None,plotOpt,'lineWid',RefreshPlot,
            size=(50,-1))
    hbox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    vbox.Add(hbox,0,wx.ALL|wx.EXPAND)

    vbox.Add((1,5))
    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add(wx.StaticText(dlg,wx.ID_ANY,' Title'),0,wx.ALL)
    plotOpt['title'] = Plot.get_title()
    w = G2G.ValidatedTxtCtrl(dlg,plotOpt,'title',OnLeave=RefreshPlot,
                                 size=(200,-1),notBlank=False)
    hbox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
    vbox.Add(hbox,0,wx.ALL|wx.EXPAND)

    for i,P in enumerate(Plot.figure.axes):
        if not P.get_visible(): continue
        if i == 0:
            lbl = ''
        else:
            lbl = str(i)
        vbox.Add((1,5))
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(wx.StaticText(dlg,wx.ID_ANY,' x label '+lbl),0,wx.ALL)
        plotOpt['xtitle'+lbl] = P.get_xaxis().get_label_text()
        w = G2G.ValidatedTxtCtrl(dlg,plotOpt,'xtitle'+lbl,OnLeave=RefreshPlot,
                                 size=(200,-1),notBlank=False)
        hbox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
        vbox.Add(hbox,0,wx.ALL|wx.EXPAND)
    
        vbox.Add((1,5))
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(wx.StaticText(dlg,wx.ID_ANY,' y label '+lbl),0,wx.ALL)
        plotOpt['ytitle'+lbl] = P.get_yaxis().get_label_text()
        w = G2G.ValidatedTxtCtrl(dlg,plotOpt,'ytitle'+lbl,OnLeave=RefreshPlot,
                                 size=(200,-1),notBlank=False)
        hbox.Add(w,0,wx.ALL|wx.ALIGN_CENTER)
        vbox.Add(hbox,0,wx.ALL|wx.EXPAND)
    
    vbox.Add((1,10),1,wx.ALL|wx.EXPAND,1)
    hbox = wx.BoxSizer(wx.HORIZONTAL)
    OKbtn = wx.Button(dlg, wx.ID_OK)
    OKbtn.Bind(wx.EVT_BUTTON,lambda event:dlg.EndModal(wx.ID_OK))
    hbox.Add((-1,-1),1,wx.ALL|wx.EXPAND,1)
    hbox.Add(OKbtn)
    hbox.Add((-1,-1),1,wx.ALL|wx.EXPAND,1)
    vbox.Add(hbox,1,wx.ALL|wx.EXPAND,1)
        
    dlg.SetSizer(vbox)
    vbox.Fit(dlg)
    #dlg.Show()
    RefreshPlot()
    dlg.ShowModal()
    
def uneqImgShow(figure,ax,Xlist,Ylist,cmap,vmin,vmax,Ylbls=[]):
    '''Plots a contour plot where point spacing varies within a dataset 
    and where the X values may differ between histograms. Note that 
    the length of Xlist and Ylist must be the same and will be the number
    of histograms to be plotted 

    :param matplotlib.figure figure:
        The figure where the plot will be placed.
    :param matplotlib.axes ax:
        The axes where the plot will be made.
    :param list Xlist:
        A list of X values for each histogram.
    :param list Ylist:
        A list of intensities for each histogram.
    :param matplotlib.colormap cmap: 
        The colormap used for shading intensities.
    :param float vmin:
        Minimum intensity.
    :param float vmax: float
        Maximum intensity.
    :param  list Ylbls: Optional.
        Label to place on each histogram. The default is [] where the axes
        are labeled normally with the first histogram numbered starting at 0.
    '''
    def midPoints(x):
        '''Return the pixel corners for a series of steps
        For the series [1,2,3,5] this will be [0.5,1.5,2.5,4,6]
        Note that n+1 points are returned for input of n points
        '''
        return np.concatenate( [[1.5*x[0] - x[1]/2], (x[:-1]+x[1:])/2, [1.5*x[-1] - x[-2]/2]] )

    lenX = len(Xlist) 
    if lenX != len(Ylist): 
        raise Exception("uneqImgShow error: unequal list lengths")
    figure.subplots_adjust(right=.85)
    #print('vmin,vmax',vmin,vmax)
    meshlist = []
    for i,(X,Y) in enumerate(zip(Xlist,Ylist)):
        #print(i,'X',min(X),max(X),'Y',min(Y),max(Y))
        meshlist.append(
            ax.pcolormesh(midPoints(X), [i-0.5,i+0.5], Y[np.newaxis,:],
                      cmap=cmap,vmin=vmin,vmax=vmax))
    # label y axis with provided labels
    if lenX == len(Ylbls):
        pos =  np.arange(lenX)
        ax.set_yticks(pos,Ylbls)
    # add the colorbar
    ax1 = figure.add_axes([0.87, 0.1, 0.04, 0.8])
    mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=mpl.colors.Normalize(vmin,vmax))
    # does not plot grid lines at present
    # if mpl.rcParams['axes.grid'] 

def initPartialOpts(phaseColors):
    '''Make sure that the options for display of partials are all defined
    '''
    for p in phaseColors:
        partialOpts[p] = partialOpts.get(p,{})
        partialOpts[p]['Show'] = partialOpts[p].get('Show',True)
        partialOpts[p]['color'] = partialOpts[p].get('color',phaseColors[p])
        partialOpts[p]['width'] = partialOpts[p].get('width','1')
        partialOpts[p]['style'] = partialOpts[p].get('style',2) # index in ltypeChoices/ltypeMPLname (see below)
        partialOpts[p]['MPLstyle'] = partialOpts[p].get('MPLstyle','--')
        
def configPartialDisplay(G2frame,phaseColors,RefreshPlot):
    '''Select which phase are diplayed with phase partials and how they are 
    displayed
    '''
    def OnSelectColor(event):
        '''Respond to a change in color
        '''
        p = event.GetEventObject().phase
        c = event.GetValue()
        # convert wx.Colour to mpl color string
        partialOpts[p]['color'] = mpl.colors.to_hex(np.array(c.Get())/255)
        phaseColors[p] = partialOpts[p]['color']
        RefreshPlot()
        
    def StyleChange(*args):
        'define MPL line style from line style index'
        for p in partialOpts:
            try: 
                partialOpts[p]['MPLstyle'] = ltypeMPLname[partialOpts[p]['style']]
            except:
                partialOpts[p]['MPLstyle'] = '-'
        RefreshPlot()

    import  wx.lib.colourselect as csel
    dlg = wx.Dialog(G2frame,
                    style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(
        wx.StaticText(dlg,wx.ID_ANY,'Set options for display of Phase partials\n(profiles for each phase)',
                          style=wx.ALIGN_CENTER),
        0,wx.ALIGN_CENTER)
    mainSizer.Add((-1,10))
    headers = ['phase name  ',' Show ',' Color ','Width','Line type']
    gsizer = wx.FlexGridSizer(cols=len(headers),hgap=2,vgap=2)
    for h in headers:
        txt = wx.StaticText(dlg,wx.ID_ANY,h)
        txt.Wrap(150)
        gsizer.Add(txt,0,wx.ALIGN_CENTER|wx.BOTTOM,6)
    initPartialOpts(phaseColors)
    for p in phaseColors:
        txt = wx.StaticText(dlg,wx.ID_ANY,p)
        txt.Wrap(150)
        gsizer.Add(txt,0,wx.ALIGN_LEFT)
        #
        ch = G2G.G2CheckBox(dlg,'',partialOpts[p],'Show',RefreshPlot)
        gsizer.Add(ch,0,wx.ALIGN_CENTER)
        #
        c = wx.Colour(mpl.colors.to_hex(partialOpts[p]['color']))
        b = csel.ColourSelect(dlg, -1, '', c)
        b.phase = p
        b.Bind(csel.EVT_COLOURSELECT, OnSelectColor)
        gsizer.Add(b,0,wx.ALL|wx.ALIGN_CENTER,3)
        #
        lwidChoices = ('0.5','0.7','1','1.5','2','2.5','3','4')
        ch = G2G.G2ChoiceButton(dlg,lwidChoices,None,None,
                                partialOpts[p],'width',RefreshPlot,
                                size=(50,-1))
        gsizer.Add(ch,0,wx.ALIGN_CENTER)
        #
        ltypeChoices = ('solid','dotted','dashed','dash-dot','dense dashed','dense dashdotted')
        ltypeMPLname = ('-',    ':',     '--',    '-.',      (0, (5, 1)),    (0, (3, 1, 1, 1)))
        ch = G2G.G2ChoiceButton(dlg,ltypeChoices,
                                partialOpts[p],'style',
                                None,None,StyleChange)
        gsizer.Add(ch,0,wx.ALIGN_CENTER)
    mainSizer.Add(gsizer)
    # OK/Cancel buttons
    btnsizer = wx.StdDialogButtonSizer()
    OKbtn = wx.Button(dlg, wx.ID_OK)
    OKbtn.SetDefault()
    btnsizer.AddButton(OKbtn)
    btnsizer.Realize()
    mainSizer.Add(btnsizer,5,wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER,0)
    mainSizer.Layout()
    dlg.SetSizer(mainSizer)
    mainSizer.Fit(dlg)
    dlg.ShowModal()
    StyleChange()
