# -*- coding: utf-8 -*-
'''Classes in :mod:`~GSASII.exports.G2export_JSON` follow:

This code is to honor my friend Robert Papoular, who wants to see what is 
inside a .gpx file.
'''
from __future__ import division, print_function
import json
import numpy as np
from .. import GSASIIfiles as G2fil

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
        elif 'G2VarObj' in str(type(obj)):
            return {"_GSASIIobj.G2VarObj":obj.varname()}
        else:
            print('Tell Brian to fix JsonEncoder to handle type=',type(obj),
                      '. Skipping for now')
            #breakpoint()
            return "sorry, I don't know how to show a {} object".format(str(type(obj)))
    
class ExportJSON(G2fil.ExportBaseclass):
    '''Implement JSON export of entire projects
    '''
    def __init__(self,G2frame):
        G2fil.ExportBaseclass.__init__(self,
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
        self.Write('[\n')
        first = True
        try:
            import wx
            wx.BeginBusyCursor()
        except:
            pass
        G2frame = self.G2frame
        # crawl through the tree, dumping as we go
        try:
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                data = []
                name = G2frame.GPXtree.GetItemText(item)
                #print('level 0',name)
                data = {name:G2frame.GPXtree.GetItemPyData(item)}
                if first:
                    first = False
                    self.Write('\n')
                else:
                    self.Write('\n, ')
                self.Write(json.dumps(
                    "=========== '{}' Tree Item ==============".format(name))+',')
                self.Write(json.dumps(data, indent=2, cls=JsonEncoder))
                item2, cookie2 = G2frame.GPXtree.GetFirstChild(item)
                while item2:
                    name2 = G2frame.GPXtree.GetItemText(item2)
                    #print('  level 1',name2)
                    self.Write(',\n')
                    self.Write(json.dumps([
                        "=========== '{}' SubItem of Tree '{}' ==============".format(name2,name)]))
                    self.Write(', ')
                    data = {name:{name2:G2frame.GPXtree.GetItemPyData(item2)}}
                    self.Write(json.dumps(data, indent=2, cls=JsonEncoder))
                    item2, cookie2 = G2frame.GPXtree.GetNextChild(item, cookie2)                            
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)                            
        finally:
            try:
                import wx
                wx.EndBusyCursor()
            except:
                pass
        self.Write(']\n')
        self.CloseFile()
