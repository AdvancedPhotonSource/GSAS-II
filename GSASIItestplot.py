# -*- coding: utf-8 -*-
#GSASIItestplot.py
'''
*GSASIItestplot: Plotting for testDeriv*
========================================

Plotting module used for script testDeriv.
'''
import wx
import wx.aui
import matplotlib as mpl
try:
    from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
except ImportError:
    from matplotlib.backends.backend_wx import FigureCanvas as Canvas
try:
    from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as Toolbar
except ImportError:
    from matplotlib.backends.backend_wxagg import Toolbar as Toolbar # name changes in wx4.0.1

class Plot(wx.Panel):
    'Creates a plotting window'
    def __init__(self, parent, id = -1, dpi = None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        self.figure = mpl.figure.Figure(dpi=dpi, #figsize=(5,7)
                                        )
        self.canvas = Canvas(self, -1, self.figure)
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar, 0 , wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)

class PlotNotebook(wx.Panel):
    'creates a Wx application and a plotting notebook'
    def __init__(self, id = -1):
        self.app = wx.App()
        self.frame = wx.Frame(None,-1,'Plotter', size=wx.Size(600,600),
            style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
        self.status = self.frame.CreateStatusBar()
        self.status.SetStatusText('Use K-box to set plot controls')
        wx.Panel.__init__(self, self.frame, id=id)
        self.nb = wx.aui.AuiNotebook(self,
            style=wx.aui.AUI_NB_DEFAULT_STYLE ^ wx.aui.AUI_NB_CLOSE_ON_ACTIVE_TAB)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)

    def Show(self):
        self.frame.Show()

    def StartEventLoop(self):
        self.Show()
        self.app.MainLoop()

    def add(self,name="plot"):
        
        def OnMotion(event):
            xpos = event.xdata
            if xpos:                                        #avoid out of frame mouse position
                ypos = event.ydata
                self.status.SetStatusText('X= %.3f Y= %.3f'%(xpos,ypos))
                
        page = Plot(self.nb)
        page.canvas.mpl_connect('motion_notify_event', OnMotion)
        self.nb.AddPage(page,name)
        return page.figure
