#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2export_JSON: ASCII .gpx Export*
------------------------------------------------------

This implements a simple exporter :class:`ExportCIF` that can export the 
contents of an entire project as a sort-of human readable (JSON) ASCII file.
This provides a way to see the contents of a GSAS-II project file, but 
this format does not provide a mechanism to change the contents of a .gpx file,
as the likelihood of breaking a data structure is too high, so these
JSON files cannot be imported back into GSAS-II. 
If you want to change the contents of a .gpx file, use :mod:`GSASIIscriptable` 
where you can access the native Python data structures and change things 
with a good chance of getting things to work. 

This code is dedicated to my friend Robert Papoular who wants to see what is 
inside a .gpx file.
'''
from __future__ import division, print_function
import json
import wx
import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIIO as G2IO

class JsonEncoder(json.JSONEncoder):
    '''This provides the ability to turn np arrays and masked arrays
    into something that json can handle.
    '''
    def default(self, obj):
        for t in [list,tuple,dict]:
            if isinstance(obj, t):
                return json.JSONEncoder.default(self, obj)
        if isinstance(obj, np.ma.core.MaskedArray):
            return {"_npmask":obj.mask.tolist(),'_npvalues':obj.tolist()}
        elif isinstance(obj, np.ndarray):
            return {"_nparray":obj.tolist()}
        else:
            print('Tell Brian to fix JsonEncoder to handle type=',type(obj),
                      '. Skipping for now')
            return "sorry, I don't know how to show a {} object".format(str(type(obj)))
    
class ExportJSON(G2IO.ExportBaseclass):
    '''Implement JSON export of entire projects
    '''
    def __init__(self,G2frame):
        G2IO.ExportBaseclass.__init__(self,
            G2frame=G2frame,
            formatName = 'ASCII JSON dump',
            extension='.json',
            longFormatName = 'Export project in ASCII a JSON dump'
            )
        self.exporttype = ['project']
        
    def Exporter(self,event=None):
        # set up for export
        self.InitExport(event)
        if self.ExportSelect(): return # set export parameters; get file name
        self.OpenFile()
        wx.BeginBusyCursor()
        G2frame = self.G2frame
        # crawl through the tree, dumping as we go
        try:
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                data = []
                name = G2frame.GPXtree.GetItemText(item)
                #print('level 0',name)
                data = {name:G2frame.GPXtree.GetItemPyData(item)}
                self.Write('\n')
                self.Write(json.dumps(
                    "=========== '{}' Tree Item ==============".format(name)))
                self.Write(json.dumps(data, indent=2, cls=JsonEncoder))
                item2, cookie2 = G2frame.GPXtree.GetFirstChild(item)
                while item2:
                    name2 = G2frame.GPXtree.GetItemText(item2)
                    #print('  level 1',name2)
                    self.Write('\n')
                    self.Write(json.dumps([
                        "=========== '{}' SubItem of Tree '{}' ==============".format(name2,name)]))
                    data = {name:{name2:G2frame.GPXtree.GetItemPyData(item)}}
                    self.Write(json.dumps(data, indent=2, cls=JsonEncoder))
                    item2, cookie2 = G2frame.GPXtree.GetNextChild(item, cookie2)                            
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)                            
        finally:
            wx.EndBusyCursor()
        self.CloseFile()
