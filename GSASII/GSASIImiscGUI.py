# -*- coding: utf-8 -*-
'''
Misc routines for GUI-based input and output, including image reading follow.

This module contains quite a bit of older code that could use some attention,
or possibly movement into other modules. It was previously called GSASIIIO.py
which is why most modules reference it as G2IO.

'''

from __future__ import division, print_function

# # Allow this to be imported without wx present.
# try:
#     import wx
# except ImportError:
#     print('wx failed')
#     # was needed by sphinx, but probably not anymore
#     class Placeholder(object):
#         def __init__(self):
#             self.Dialog = object
#     wx = Placeholder()
import math
import os
import re
import copy
import platform
import pickle
import sys
import random as ran

import numpy as np
import numpy.ma as ma
import wx

from . import GSASIIpath
from . import GSASIIdataGUI as G2gd
from . import GSASIIobj as G2obj
#import GSASIIpwdGUI as G2pdG
from . import GSASIIimgGUI as G2imG
from . import GSASIIElem as G2el
from . import GSASIIfiles as G2fil
from . import GSASIIctrlGUI as G2G
from . import GSASIImath as G2mth
from . import GSASIIElem as G2elem
from . import GSASIIspc as G2spc
from . import GSASIIlattice as G2lat
from . import GSASIIpwd as G2pwd

DEBUG = False       #=True for various prints
TRANSP = False      #=true to transpose images for testing
if GSASIIpath.GetConfigValue('Transpose'): TRANSP = True
npsind = lambda x: np.sin(x*np.pi/180.)

def FileDlgFixExt(dlg,file):
    'this is needed to fix a problem in linux wx.FileDialog'
    ext = dlg.GetWildcard().split('|')[2*dlg.GetFilterIndex()+1].strip('*')
    if ext not in file:
        file += ext
    return file

def GetPowderPeaks(fileName):
    'Read powder peaks from a file'
    sind = lambda x: math.sin(x*math.pi/180.)
    asind = lambda x: 180.*math.asin(x)/math.pi
    wave = 1.54052
    File = open(fileName,'r')
    Comments = []
    peaks = []
    S = File.readline()
    while S:
        if S[:1] == '#':
            Comments.append(S[:-1])
        else:
            item = S.split()
            if len(item) == 1:
                peaks.append([float(item[0]),1.0])
            elif len(item) > 1:
                peaks.append([float(item[0]),float(item[1])])
        S = File.readline()
    File.close()
    if Comments:
       print ('Comments on file:')
       for Comment in Comments:
            print (Comment)
            if 'wavelength' in Comment:
                wave = float(Comment.split('=')[1])
    Peaks = []
    if peaks[0][0] > peaks[-1][0]:          # d-spacings - assume CuKa
        for peak in peaks:
            dsp = peak[0]
            sth = wave/(2.0*dsp)
            if sth < 1.0:
                tth = 2.0*asind(sth)
            else:
                break
            Peaks.append([tth,peak[1],True,False,0,0,0,dsp,0.0])
    else:                                   #2-thetas - assume Cuka (for now)
        for peak in peaks:
            tth = peak[0]
            dsp = wave/(2.0*sind(tth/2.0))
            Peaks.append([tth,peak[1],True,False,0,0,0,dsp,0.0])
    limits = [1000.,0.]
    for peak in Peaks:
        limits[0] = min(limits[0],peak[0])
        limits[1] = max(limits[1],peak[0])
    limits[0] = max(1.,(int(limits[0]-1.)/5)*5.)
    limits[1] = min(170.,(int(limits[1]+1.)/5)*5.)
    return Comments,Peaks,limits,wave

def GetCheckImageFile(G2frame,treeId):
    '''Try to locate an image file if the project and image have been moved
    together. If the image file cannot be found, request the location from
    the user.

    :param wx.Frame G2frame: main GSAS-II Frame and data object
    :param wx.Id treeId: Id for the main tree item for the image
    :returns: Npix,imagefile,imagetag with (Npix) number of pixels,
       imagefile, if it exists, or the name of a file that does exist or False if the user presses Cancel
       and (imagetag) an optional image number

    '''
    Npix,Imagefile,imagetag = G2frame.GPXtree.GetImageLoc(treeId)
    if isinstance(Imagefile,list):
        imagefile,imagetag = Imagefile
    else:
        imagefile = Imagefile
    if not os.path.exists(imagefile):
        print ('Image file '+imagefile+' not found')
        fil = imagefile.replace('\\','/') # windows?!
        # see if we can find a file with same name or in a similarly named sub-dir
        pth,fil = os.path.split(fil)
        prevpth = None
        while pth and pth != prevpth:
            prevpth = pth
            if os.path.exists(os.path.join(G2frame.dirname,fil)):
                print ('found image file '+os.path.join(G2frame.dirname,fil))
                imagefile = os.path.join(G2frame.dirname,fil)
                G2frame.GPXtree.UpdateImageLoc(treeId,imagefile)
                return Npix,imagefile,imagetag
            pth,enddir = os.path.split(pth)
            fil = os.path.join(enddir,fil)
        # not found as a subdirectory, drop common parts of path for last saved & image file names
        #    if image was .../A/B/C/imgs/ima.ge
        #      & GPX was  .../A/B/C/refs/fil.gpx but is now .../NEW/TEST/TEST1
        #    will look for .../NEW/TEST/TEST1/imgs/ima.ge, .../NEW/TEST/imgs/ima.ge, .../NEW/imgs/ima.ge and so on
        Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        gpxPath = Controls.get('LastSavedAs','').replace('\\','/').split('/') # blank in older .GPX files
        imgPath = imagefile.replace('\\','/').split('/')
        for p1,p2 in zip(gpxPath,imgPath):
            if p1 == p2:
                gpxPath.pop(0),imgPath.pop(0)
            else:
                break
        fil = os.path.join(*imgPath) # file with non-common prefix elements
        prevpth = None
        pth = os.path.abspath(G2frame.dirname)
        while pth and pth != prevpth:
            prevpth = pth
            if os.path.exists(os.path.join(pth,fil)):
                print ('found image file '+os.path.join(pth,fil))
                imagefile = os.path.join(pth,fil)
                G2frame.GPXtree.UpdateImageLoc(treeId,imagefile)
                return Npix,imagefile,imagetag
            pth,enddir = os.path.split(pth)
        #GSASIIpath.IPyBreak()

    if not os.path.exists(imagefile):
        # note that this fails (at least on Mac) to get an image during the GUI initialization
        prevnam = os.path.split(imagefile)[1]
        prevext = os.path.splitext(imagefile)[1]
        wildcard = 'Image format (*'+prevext+')|*'+prevext
        dlg = wx.FileDialog(G2frame, 'Previous image file ('+prevnam+') not found; open here', '.', prevnam,
                            wildcard,wx.FD_OPEN)
        try:
            dlg.SetFilename(''+os.path.split(imagefile)[1])
            if dlg.ShowModal() == wx.ID_OK:
                imagefile = dlg.GetPath()
                G2frame.GPXtree.UpdateImageLoc(treeId,imagefile)
            else:
                imagefile = None # was False
        finally:
            dlg.Destroy()
    return Npix,imagefile,imagetag

def EditImageParms(parent,Data,Comments,Image,filename):
    dlg = wx.Dialog(parent, wx.ID_ANY, 'Edit image parameters',
                    style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
    def onClose(event):
        dlg.EndModal(wx.ID_OK)
    mainsizer = wx.BoxSizer(wx.VERTICAL)
    h,w = Image.size[:2]
    mainsizer.Add(wx.StaticText(dlg,wx.ID_ANY,'File '+str(filename)+'\nImage size: '+str(h)+' x '+str(w)),
        0,wx.ALIGN_LEFT|wx.ALL, 2)

    vsizer = wx.BoxSizer(wx.HORIZONTAL)
    vsizer.Add(wx.StaticText(dlg,wx.ID_ANY,u'Wavelength (\xC5) '),
        0,wx.ALIGN_LEFT|wx.ALL, 2)
    wdgt = G2G.ValidatedTxtCtrl(dlg,Data,'wavelength')
    vsizer.Add(wdgt)
    mainsizer.Add(vsizer,0,wx.ALIGN_LEFT|wx.ALL, 2)

    vsizer = wx.BoxSizer(wx.HORIZONTAL)
    vsizer.Add(wx.StaticText(dlg,wx.ID_ANY,u'Pixel size (\xb5m). Width '),
        0,wx.ALIGN_LEFT|wx.ALL, 2)
    wdgt = G2G.ValidatedTxtCtrl(dlg,Data['pixelSize'],0,size=(50,-1))
    vsizer.Add(wdgt)
    vsizer.Add(wx.StaticText(dlg,wx.ID_ANY,u'  Height '),wx.ALIGN_LEFT|wx.ALL, 2)
    wdgt = G2G.ValidatedTxtCtrl(dlg,Data['pixelSize'],1,size=(50,-1))
    vsizer.Add(wdgt)
    mainsizer.Add(vsizer,0,wx.ALIGN_LEFT|wx.ALL, 2)

    vsizer = wx.BoxSizer(wx.HORIZONTAL)
    vsizer.Add(wx.StaticText(dlg,wx.ID_ANY,u'Sample to detector (mm) '),
        0,wx.ALIGN_LEFT|wx.ALL, 2)
    wdgt = G2G.ValidatedTxtCtrl(dlg,Data,'distance')
    vsizer.Add(wdgt)
    mainsizer.Add(vsizer,0,wx.ALIGN_LEFT|wx.ALL, 2)

    vsizer = wx.BoxSizer(wx.HORIZONTAL)
    vsizer.Add(wx.StaticText(dlg,wx.ID_ANY,u'Beam center (pixels). X = '),
        0,wx.ALIGN_LEFT|wx.ALL, 2)
    wdgt = G2G.ValidatedTxtCtrl(dlg,Data['center'],0,size=(75,-1))
    vsizer.Add(wdgt)
    vsizer.Add(wx.StaticText(dlg,wx.ID_ANY,u'  Y = '),wx.ALIGN_LEFT|wx.ALL, 2)
    wdgt = G2G.ValidatedTxtCtrl(dlg,Data['center'],1,size=(75,-1))
    vsizer.Add(wdgt)
    mainsizer.Add(vsizer,0,wx.ALIGN_LEFT|wx.ALL, 2)

    vsizer = wx.BoxSizer(wx.HORIZONTAL)
    vsizer.Add(wx.StaticText(dlg,wx.ID_ANY,u'Comments '),
        0,wx.ALIGN_LEFT|wx.ALL, 2)
    wdgt = G2G.ValidatedTxtCtrl(dlg,Comments,0,size=(250,-1))
    vsizer.Add(wdgt)
    mainsizer.Add(vsizer,0,wx.ALIGN_LEFT|wx.ALL, 2)

    btnsizer = wx.StdDialogButtonSizer()
    OKbtn = wx.Button(dlg, wx.ID_OK, 'Continue')
    OKbtn.SetDefault()
    OKbtn.Bind(wx.EVT_BUTTON,onClose)
    btnsizer.AddButton(OKbtn) # not sure why this is needed
    btnsizer.Realize()
    mainsizer.Add(btnsizer, 1, wx.ALL|wx.EXPAND, 5)
    dlg.SetSizer(mainsizer)
    dlg.CenterOnParent()
    dlg.ShowModal()

def LoadImage2Tree(imagefile,G2frame,Comments,Data,Npix,Image):
    '''Load an image into the tree. Saves the location of the image, as well as the
    ImageTag (where there is more than one image in the file), if defined.
    '''
    ImgNames = G2gd.GetGPXtreeDataNames(G2frame,['IMG ',])
    TreeLbl = 'IMG '+os.path.basename(imagefile)
    ImageTag = Data.get('ImageTag')
    ImageSection = Data.get('ImageSection','')  #used only in HDF5 at present
    if ImageTag:
        if ImageSection and ImageTag[1] is None:
            TreeLbl += f" {ImageSection}"
        elif ImageSection:
            TreeLbl += f" {ImageSection} #{ImageTag[1]:04d}"
        else:
            TreeLbl += f' #{ImageTag:04d}'
        imageInfo = (imagefile,ImageTag)
    else:
        imageInfo = imagefile
    TreeName = G2obj.MakeUniqueLabel(TreeLbl,ImgNames)
    Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=TreeName)
    G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Comments'),Comments)
    Imax = np.amax(Image)
    if G2frame.imageDefault:
        Data.update(copy.deepcopy(G2frame.imageDefault))
        Data['showLines'] = True
        Data['ring'] = []
        Data['rings'] = []
        Data['cutoff'] = 10.
        Data['pixLimit'] = 20
        Data['edgemin'] = 100000000
        Data['calibdmin'] = 0.5
        Data['calibskip'] = 0
        Data['ellipses'] = []
        Data['calibrant'] = ''
        Data['GonioAngles'] = [0.,0.,0.]
        Data['DetDepthRef'] = False
    else:
        Data['type'] = 'PWDR'
        Data['color'] = GSASIIpath.GetConfigValue('Contour_color','Paired')
        if 'tilt' not in Data:          #defaults if not preset in e.g. Bruker importer
            Data['tilt'] = 0.0
            Data['rotation'] = 0.0
            Data['pixLimit'] = 20
            Data['calibdmin'] = 0.5
            Data['cutoff'] = 10.
        Data['showLines'] = False
        Data['calibskip'] = 0
        Data['ring'] = []
        Data['rings'] = []
        Data['edgemin'] = 100000000
        Data['ellipses'] = []
        Data['GonioAngles'] = [0.,0.,0.]
        Data['DetDepth'] = 0.
        Data['DetDepthRef'] = False
        Data['calibrant'] = ''
        Data['IOtth'] = [5.0,50.0]
        if GSASIIpath.GetConfigValue('Image_2theta_min'):
            try:
                Data['IOtth'][0] = float(GSASIIpath.GetConfigValue('Image_2theta_min'))
            except:
                pass
        if GSASIIpath.GetConfigValue('Image_2theta_max'):
            try:
                Data['IOtth'][1] = float(GSASIIpath.GetConfigValue('Image_2theta_max'))
            except:
                pass
        Data['LRazimuth'] = [0.,180.]
        Data['azmthOff'] = 0.0
        Data['outChannels'] = 2500
        Data['outAzimuths'] = 1
        Data['centerAzm'] = False
        Data['fullIntegrate'] = GSASIIpath.GetConfigValue('fullIntegrate',True)
        Data['setRings'] = False
        Data['background image'] = ['',-1.0]
        Data['dark image'] = ['',-1.0]
        Data['Flat Bkg'] = 0.0
        Data['Oblique'] = [0.5,False]
    Data['setDefault'] = False
    Data['range'] = [(0,Imax),[0,Imax]]
    G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Image Controls'),Data)
    Masks = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[],
                 'Thresholds':[(0,Imax),[0,Imax]],
                 'SpotMask':{'esdMul':3.,'spotMask':None}}
    G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Masks'),Masks)
    G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Stress/Strain'),
        {'Type':'True','d-zero':[],'Sample phi':0.0,'Sample z':0.0,'Sample load':0.0})
    G2frame.GPXtree.SetItemPyData(Id,[Npix,imageInfo])
    G2frame.PickId = Id
    G2frame.PickIdText = G2frame.GetTreeItemsList(G2frame.PickId)
    G2frame.Image = Id

def ReadImages(G2frame,imagefile):
    '''Read one or more images from a file and put them into the Tree
    using image importers. Called only in :meth:`AutoIntFrame.OnTimerLoop`.

    ToDo: Images are most commonly read in :meth:`GSASIIdataGUI.GSASII.OnImportGeneric`
    which is called from :meth:`GSASIIdataGUI.GSASII.OnImportImage`
    it would be good if these routines used a common code core so that changes need to
    be made in only one place.

    :param wx.Frame G2frame: main GSAS-II Frame and data object.
    :param str imagefile: name of image file

    :returns: a list of the id's of the IMG tree items created
    '''
    # determine which formats are compatible with this file
    primaryReaders = []
    secondaryReaders = []
    for rd in G2frame.ImportImageReaderlist:
        flag = rd.ExtensionValidator(imagefile)
        if flag is None:
            secondaryReaders.append(rd)
        elif flag:
            primaryReaders.append(rd)
    if len(secondaryReaders) + len(primaryReaders) == 0:
        print('Error: No matching format for file '+imagefile)
        raise Exception('No image read')
    errorReport = ''
    rdbuffer = {} # create temporary storage for file reader
    for rd in primaryReaders+secondaryReaders:
        rd.ReInitialize() # purge anything from a previous read
        rd.errors = "" # clear out any old errors
        if not rd.ContentsValidator(imagefile): # rejected on cursory check
            errorReport += "\n  "+rd.formatName + ' validator error'
            if rd.errors:
                errorReport += ': '+rd.errors
                continue
        ParentFrame = G2frame
        block = 0
        repeat = True
        CreatedIMGitems = []
        while repeat: # loop if the reader asks for another pass on the file
            block += 1
            repeat = False
            if GSASIIpath.GetConfigValue('debug'):
                flag = rd.Reader(imagefile,ParentFrame,blocknum=block,Buffer=rdbuffer)
            else:
                flag = False
                try:
                    flag = rd.Reader(imagefile,ParentFrame,blocknum=block,Buffer=rdbuffer)
                except rd.ImportException as detail:
                    rd.errors += "\n  Read exception: "+str(detail)
                except Exception as detail:
                    import traceback
                    rd.errors += "\n  Unhandled read exception: "+str(detail)
                    rd.errors += "\n  Traceback info:\n"+str(traceback.format_exc())
            if flag: # this read succeeded
                if rd.Image is None:
                    raise Exception('No image read. Strange!')
                if GSASIIpath.GetConfigValue('Transpose'):
                    print ('Transposing Image!')
                    rd.Image = rd.Image.T
                rd.Data['ImageTag'] = rd.repeatcount
                rd.readfilename = imagefile
                # Load generic metadata, as configured
                G2fil.GetColumnMetadata(rd)
                LoadImage2Tree(imagefile,G2frame,rd.Comments,rd.Data,rd.Npix,rd.Image)
                repeat = rd.repeat
            CreatedIMGitems.append(G2frame.Image)
        if CreatedIMGitems: return CreatedIMGitems
    else:
        print('Error reading file '+imagefile)
        print('Error messages(s)\n'+errorReport)
        return []
        #raise Exception('No image read')

def SaveMultipleImg(G2frame):
    if not G2frame.GPXtree.GetCount():
        print ('no images!')
        return
    choices = G2gd.GetGPXtreeDataNames(G2frame,['IMG ',])
    if len(choices) == 1:
        names = choices
    else:
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Stress/Strain fitting','Select images to fit:',choices)
        dlg.SetSelections([])
        names = []
        if dlg.ShowModal() == wx.ID_OK:
            names = [choices[sel] for sel in dlg.GetSelections()]
        dlg.Destroy()
    if not names: return
    for name in names:
        Id = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, name)
        Npix,imagefile,imagetag = G2frame.GPXtree.GetImageLoc(Id)
        imroot = os.path.splitext(imagefile)[0]
        if imagetag:
            imroot += '_' + str(imagetag)
        Data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'Image Controls'))
        print('Writing '+imroot+'.imctrl')
        File = open(imroot+'.imctrl','w')
        keys = ['type','wavelength','calibrant','distance','center',
                    'tilt','rotation','azmthOff','fullIntegrate','LRazimuth',
                    'IOtth','outChannels','outAzimuths','invert_x','invert_y','DetDepth',
                    'calibskip','pixLimit','cutoff','calibdmin','chisq','Flat Bkg',
                    'binType','SampleShape','PolaVal','SampleAbs','dark image','background image']
        for key in keys:
            if key not in Data: continue    #uncalibrated!
            File.write(key+':'+str(Data[key])+'\n')
        File.close()
        mask = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'Masks'))
        G2imG.CleanupMasks(mask)
        print('Writing '+imroot+'.immask')
        File = open(imroot+'.immask','w')
        for key in ['Points','Rings','Arcs','Polygons','Frames','Thresholds']:
            File.write(key+':'+str(mask[key])+'\n')
        File.close()

def PutG2Image(filename,Comments,Data,Npix,image):
    'Write an image as a python pickle - might be better as an .edf file?'
    File = open(filename,'wb')
    pickle.dump([Comments,Data,Npix,image],File,2)
    File.close()
    return

objectScanIgnore = [int,bool,float,str,np.float64,np.float32,np.int32,np.int64,np.int16,np.ndarray,G2obj.G2VarObj,G2obj.ExpressionObj,np.bool_]
try:
    objectScanIgnore += [ma.MaskedArray] # fails in doc builds
except AttributeError:
    pass

def objectScan(data,tag,indexStack=[]):
    '''Recursively scan an object looking for unexpected data types.
    This is used in debug mode to scan .gpx files for objects we did not
    intend to be there.
    '''
    if type(data) is list or type(data) is tuple:
        for i in range(len(data)):
            val = objectScan(data[i],tag,indexStack+[i])
            if val:
                data[i] = val
                print('...fixed')
    elif type(data) is dict:
        for key in data:
            val = objectScan(data[key],tag,indexStack+[key])
            if val:
                data[key] = val
                print('...fixed')
    elif data is None:
        return None
    elif type(data) in objectScanIgnore:
        return None
    # not always recognized:
    elif 'GSASIIobj.G2VarObj' in str(type(data)):
        return None
    elif 'GSASIIobj.ExpressionObj' in str(type(data)):
        return None
    else:
        s = 'unexpected object in '+tag
        for i in indexStack:
            s += "[{}]".format(i)
        #print(s,data.__class__.__name__) # loses full name of class
        print(s,type(data))
        global unexpectedObject
        unexpectedObject = True
        # fix bad objects
        if "gdi.Colour" in str(type(data)):
            return tuple(data)
        return

def pickleLoad(fp):
    return pickle.load(fp,encoding='latin-1')

def ProjFileOpen(G2frame,showProvenance=True):
    'Read a GSAS-II project file and load into the G2 data tree'
    if not os.path.exists(G2frame.GSASprojectfile):
        print ('\n*** Error attempt to open project file that does not exist:\n   '+
               str(G2frame.GSASprojectfile))
        return
    LastSavedUsing = None
    filep = open(G2frame.GSASprojectfile,'rb')
    if showProvenance: print ('loading from file: '+G2frame.GSASprojectfile)
    GPXphase = os.path.splitext(G2frame.GSASprojectfile)[0]+'.seqPhase'
    GPXhist = os.path.splitext(G2frame.GSASprojectfile)[0]+'.seqHist'
    deleteSeq = False
    hist = None
    tmpHistIndex = {}
    updateFromSeq = False
    if os.path.exists(GPXphase) and os.path.exists(GPXhist):
        dlg = wx.MessageDialog(G2frame,
            'Load results from crashed sequential fit?\nNo deletes the files!', 'Recover partial sequential fit?', wx.YES | wx.NO | wx.CANCEL)
        dlg.CenterOnParent()
        try:
            result = dlg.ShowModal()
            deleteSeq = result != wx.ID_CANCEL
            if result == wx.ID_YES:
                updateFromSeq = True
                fp = open(GPXphase,'rb')
                data = pickleLoad(fp) # first block in file should be Phases
                if data[0][0] != 'Phases':
                    raise Exception('Unexpected block in {} file. How did this happen?'
                            .format(GPXphase))
                Phases = {}
                for name,vals in data[1:]:
                    Phases[name] = vals
                name,CovData = pickleLoad(fp)[0] # 2nd block in file should be Covariance
                name,RigidBodies = pickleLoad(fp)[0] # 3rd block in file should be Rigid Bodies
                fp.close()
                # index the histogram updates
                hist = open(GPXhist,'rb')
                try:
                    while True:
                        loc = hist.tell()
                        datum = pickleLoad(hist)[0]
                        tmpHistIndex[datum[0]] = loc
                except EOFError:
                    pass
        finally:
            dlg.Destroy()
    wx.BeginBusyCursor()
    try:
        if GSASIIpath.GetConfigValue('show_gpxSize'):
            posPrev = 0
            sizeList = {}
        while True:
            try:
                data = pickleLoad(filep)
            except EOFError:
                break
            datum = data[0]
            if GSASIIpath.GetConfigValue('show_gpxSize'):
                sizeList[datum[0]] = filep.tell()-posPrev
                posPrev = filep.tell()
            # scan the GPX file for unexpected objects
            if GSASIIpath.GetConfigValue('debug'):
                global unexpectedObject
                unexpectedObject = False
                objectScan(data,'tree item "{}" entry '.format(datum[0]))
                #if unexpectedObject:
                #    print(datum[0])
                #    GSASIIpath.IPyBreak()
            Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=datum[0])
            if datum[0] == 'Phases' and GSASIIpath.GetConfigValue('SeparateHistPhaseTreeItem',False):
                G2frame.GPXtree.AppendItem(parent=G2frame.root,text='Hist/Phase')
            if updateFromSeq and datum[0] == 'Phases':
                for pdata in data[1:]:
                    if pdata[0] in Phases:
                        pdata[1].update(Phases[pdata[0]])
            elif updateFromSeq and datum[0] == 'Covariance':
                data[0][1] = CovData
            elif updateFromSeq and datum[0] == 'Rigid bodies':
                data[0][1] = RigidBodies
            elif updateFromSeq and datum[0] in tmpHistIndex:
                hist.seek(tmpHistIndex[datum[0]])
                hdata = pickleLoad(hist)
                if data[0][0] != hdata[0][0]:
                    print('Error! Updating {} with {}'.format(data[0][0],hdata[0][0]))
                datum = hdata[0]
                xferItems = ['Background','Instrument Parameters','Sample Parameters','Reflection Lists']
                hItems = {name:j+1 for j,(name,val) in enumerate(hdata[1:]) if name in xferItems}
                for j,(name,val) in enumerate(data[1:]):
                    if name not in xferItems: continue
                    data[j+1][1] = hdata[hItems[name]][1]
            if datum[0].startswith('PWDR'):
                if 'ranId' not in datum[1][0]: # patch: add random Id if not present
                    datum[1][0]['ranId'] = ran.randint(0,sys.maxsize)
                G2frame.GPXtree.SetItemPyData(Id,datum[1][:3])  #temp. trim off junk (patch?)
            elif datum[0].startswith('HKLF'):
                if 'ranId' not in datum[1][0]: # patch: add random Id if not present
                    datum[1][0]['ranId'] = ran.randint(0,sys.maxsize)
                G2frame.GPXtree.SetItemPyData(Id,datum[1])
            else:
                G2frame.GPXtree.SetItemPyData(Id,datum[1])
                if datum[0] == 'Controls' and 'LastSavedUsing' in datum[1]:
                    LastSavedUsing = datum[1]['LastSavedUsing']
                if datum[0] == 'Controls' and 'PythonVersions' in datum[1] and GSASIIpath.GetConfigValue('debug') and showProvenance:
                    print('DBG_Packages used to create .GPX file:')
                    if 'dict' in str(type(datum[1]['PythonVersions'])):  #patch
                        for p in sorted(datum[1]['PythonVersions'],key=lambda s: s.lower()):
                            print("  {:<14s}: {:s}".format(p[0],p[1]))
                    else:
                        for p in datum[1]['PythonVersions']:
                            print("  {:<12s} {:s}".format(p[0]+':',p[1]))
            oldPDF = False
            for datus in data[1:]:
#patch - 1/23/17 PDF cleanup
                if datus[0][:4] in ['I(Q)','S(Q)','F(Q)','G(R)']:
                    oldPDF = True
                    data[1][1][datus[0][:4]] = copy.deepcopy(datus[1][:2])
                    continue
#end PDF cleanup
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
            if 'PDF ' in datum[0][:4] and oldPDF:
                sub = G2frame.GPXtree.AppendItem(Id,'PDF Peaks')
                G2frame.GPXtree.SetItemPyData(sub,{'Limits':[1.,5.],'Background':[2,[0.,-0.2*np.pi],False],'Peaks':[]})
            if datum [0].startswith('IMG'):                   #retrieve image default flag & data if set
                Data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Image Controls'))
                if Data['setDefault']:
                    G2frame.imageDefault = Data
                    G2frame.imageDefault['setDefault'] = False
                    if 'formatName' in G2frame.imageDefault: del G2frame.imageDefault['formatName']
        if LastSavedUsing:
            print('GPX load successful. Last saved with GSAS-II revision '+LastSavedUsing)
        else:
            print('project load successful')
        if GSASIIpath.GetConfigValue('show_gpxSize'):
            print(50*'=')
            print('File section sizes (Kb)')
            for item in sizeList:
                print('  {:20s} {:10.3f}'.format(
                    item[:20],sizeList[item]/1024.))
            print(50*'=')
        G2frame.NewPlot = True
    except Exception as errmsg:
        if GSASIIpath.GetConfigValue('debug'):
            print('\nError reading GPX file:',errmsg)
            import traceback
            print (traceback.format_exc())
        msg = wx.MessageDialog(G2frame,message="Error reading file "+
            str(G2frame.GSASprojectfile)+". This is not a current GSAS-II .gpx file",
            caption="Load Error",style=wx.ICON_ERROR | wx.OK | wx.STAY_ON_TOP)
        msg.ShowModal()
    finally:
        filep.close()
        wx.EndBusyCursor()
        G2frame.Status.SetStatusText('Mouse RB drag/drop to reorder',0)
    if deleteSeq:
        if hist: hist.close()
        try:
            os.remove(GPXphase)
        except:
            print('Warning: unable to delete {}'.format(GPXphase))
        try:
            os.remove(GPXhist)
        except:
            print('Warning: unable to delete {}'.format(GPXhist))
    G2frame.SetTitleByGPX()
    if LastSavedUsing:
        try:
            G2G.updateNotifier(G2frame,int(LastSavedUsing.split()[0]))
        except:
            pass

def ProjFileSave(G2frame):
    'Save a GSAS-II project file'
    if not G2frame.GPXtree.IsEmpty():
        try:
            file = open(G2frame.GSASprojectfile,'wb')
        except PermissionError:
            G2G.G2MessageBox(G2frame,'Read only file','Project cannot be saved; change permission & try again')
            return
        print ('save to file: '+G2frame.GSASprojectfile)
        # stick the file name into the tree and version info into tree so they are saved.
        # (Controls should always have been created in tree at this point)
        try:
            Controls = G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
            Controls['PythonVersions'] = G2frame.PackageVersions
            Controls['LastSavedAs'] = os.path.abspath(G2frame.GSASprojectfile)
            Controls['LastSavedUsing'] = f"#{GSASIIpath.GetVersionNumber()}, {GSASIIpath.GetVersionTag()}"
            if GSASIIpath.HowIsG2Installed().startswith('git'):
                g2repo = GSASIIpath.openGitRepo(GSASIIpath.path2GSAS2)
                commit = g2repo.head.commit
                Controls['LastSavedUsing'] += f" git {commit.hexsha[:8]}"
            else:
                gv = getSavedVersionInfo()
                if gv is not None:
                    Controls['LastSavedUsing'] += f" static {gv.git_version[:8]}"
        except:
            pass
        wx.BeginBusyCursor()
        try:
            item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
            while item:
                data = []
                name = G2frame.GPXtree.GetItemText(item)
                if name.startswith('Hist/Phase'):  # skip over this
                    item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
                    continue
                data.append([name,G2frame.GPXtree.GetItemPyData(item)])
                item2, cookie2 = G2frame.GPXtree.GetFirstChild(item)
                while item2:
                    name = G2frame.GPXtree.GetItemText(item2)
                    data.append([name,G2frame.GPXtree.GetItemPyData(item2)])
                    item2, cookie2 = G2frame.GPXtree.GetNextChild(item, cookie2)
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
                pickle.dump(data,file,2)
            file.close()
            pth = os.path.split(os.path.abspath(G2frame.GSASprojectfile))[0]
            if GSASIIpath.GetConfigValue('Save_paths'): G2G.SaveGPXdirectory(pth)
            G2frame.LastGPXdir = pth
        finally:
            wx.EndBusyCursor()
        print('project save successful')

def SaveIntegration(G2frame,PickId,data,Overwrite=False):
    'Save image integration results as powder pattern(s)'
    waves = {'Cu':[1.54051,1.54433],'Ti':[2.74841,2.75207],'Cr':[2.28962,2.29351],
        'Fe':[1.93597,1.93991],'Co':[1.78892,1.79278],'Mo':[0.70926,0.713543],
        'Ag':[0.559363,0.563775]}
    azms = G2frame.Integrate[1]
    X = G2frame.Integrate[2][:-1]
    N = len(X)
    Id = G2frame.GPXtree.GetItemParent(PickId)
    name = G2frame.GPXtree.GetItemText(Id)
    name = name.replace('IMG ',data['type']+' ')
    Comments = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id, 'Comments'))
    Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    Comments.append('Dark image = %s\n'%str(data['dark image']))
    Comments.append('Background image = %s\n'%str(data['background image']))
    Comments.append('Gain map = %s\n'%str(data['Gain map']))

    if 'PWDR' in name:
        if 'target' in data:
            names = ['Type','Lam1','Lam2','I(L2)/I(L1)','Zero','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth']
            codes = [0 for i in range(14)]
        else:
            if data.get('IfPink',False):
                names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','Z','alpha-0','alpha-1','beta-0','beta-1','Azimuth']
                codes = [0 for i in range(15)]
            else:
                names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth']
                codes = [0 for i in range(12)]
    elif 'SASD' in name:
        names = ['Type','Lam','Zero','Azimuth']
        codes = [0 for i in range(4)]
        X = 4.*np.pi*npsind(X/2.)/data['wavelength']    #convert to q
    Xminmax = [X[0],X[-1]]
    Azms = np.zeros(data['outAzimuths'])
    dazm = 0.
    if data['outAzimuths'] > 1:
        dazm = np.min(np.abs(np.diff(azms)))/2.
    G2frame.IntgOutList = []
    for i,azm in enumerate(azms[:-1]):
        Aname = name+" Azm= %.2f"%((azm+dazm)%360.)
        item, cookie = G2frame.GPXtree.GetFirstChild(G2frame.root)
        # if Overwrite delete any duplicate
        if Overwrite and G2gd.GetGPXtreeItemId(G2frame,G2frame.root,Aname):
            print('Replacing '+Aname)
            item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,Aname)
            G2frame.GPXtree.Delete(item)
        else:
            nOcc = 0
            while item:
                Name = G2frame.GPXtree.GetItemText(item)
                if Aname in Name:
                    nOcc += 1
                item, cookie = G2frame.GPXtree.GetNextChild(G2frame.root, cookie)
            if nOcc:
                Aname += '(%d)'%(nOcc)
        Sample = G2obj.SetDefaultSample()       #set as Debye-Scherrer
        Sample['Gonio. radius'] = data['distance']
        Sample['Omega'] = data['GonioAngles'][0]
        Sample['Chi'] = data['GonioAngles'][1]
        Sample['Phi'] = data['GonioAngles'][2]
        Sample['Azimuth'] = (azm+dazm)%360.    #put here as bin center
        polariz = data['PolaVal'][0]
        for item in Comments:
            for key in ('Temperature','Pressure','Time','FreePrm1','FreePrm2','FreePrm3','Omega',
                'Chi','Phi'):
                if key.lower() in item.lower():
                    try:
                        Sample[key] = float(item.split('=')[1])
                    except:
                        pass
            if 'label_prm' in item.lower():
                for num in ('1','2','3'):
                    if 'label_prm'+num in item.lower():
                        Controls['FreePrm'+num] = item.split('=')[1].strip()
        if 'PWDR' in Aname:
            if 'target' in data:    #from lab x-ray 2D imaging data
                wave1,wave2 = waves[data['target']]
                parms = ['PXC',wave1,wave2,0.5,0.0,polariz,290.,-40.,30.,6.,-14.,0.0,0.0001,Azms[i]]
            else:
                if data.get('IfPink',False):
                    parms = ['PXB',data['wavelength'],0.0,polariz,0.,8000.,-150.,-24.,0.,0.,0.,13.,-1300.,3.,-7.,Azms[i]]   #from Sect 35 LSS
                else:
                    parms = ['PXC',data['wavelength'],0.0,polariz,1.0,-0.10,0.4,0.30,1.0,0.0,0.0001,Azms[i]]
        elif 'SASD' in Aname:
            Sample['Trans'] = data['SampleAbs'][0]
            parms = ['LXC',data['wavelength'],0.0,Azms[i]]
        Y = G2frame.Integrate[0][i]
        Ymin = np.min(Y)
        Ymax = np.max(Y)
        W = np.where(Y>0.,1./Y,1.e-6)                    #probably not true
        Id = G2frame.GPXtree.AppendItem(parent=G2frame.root,text=Aname)
        G2frame.IntgOutList.append(Id)
        G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Comments'),Comments)
        G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Limits'),copy.deepcopy([tuple(Xminmax),Xminmax]))
        if 'PWDR' in Aname:
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Background'),[['chebyschev-1',1,3,1.0,0.0,0.0],
                {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[],'background PWDR':['',1.0,False]}])
        inst = [dict(zip(names,zip(parms,parms,codes))),{}]
        for item in inst[0]:
            inst[0][item] = list(inst[0][item])
        G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Instrument Parameters'),inst)
        if 'PWDR' in Aname:
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Sample Parameters'),Sample)
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Peak List'),{'sigDict':{},'peaks':[]})
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Index Peak List'),[[],[]])
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Unit Cells List'),[])
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Reflection Lists'),{})
        elif 'SASD' in Aname:
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Substances'),G2pwd.SetDefaultSubstances())
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Sample Parameters'),Sample)
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Models'),G2pwd.SetDefaultSASDModel())
        valuesdict = {
            'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxsize),'Offset':[0.0,0.0],'delOffset':0.02*Ymax,
            'refOffset':-0.1*Ymax,'refDelt':0.1*Ymax,'Yminmax':[Ymin,Ymax]}
        G2frame.GPXtree.SetItemPyData(Id,[valuesdict,
            [np.array(X),np.array(Y),np.array(W),np.zeros(N),np.zeros(N),np.zeros(N)]])
    return Id       #last powder pattern generated

def XYsave(G2frame,XY,labelX='X',labelY='Y',names=[]):
    'Save XY table data'
    pth = G2G.GetExportPath(G2frame)
    dlg = wx.FileDialog(
        G2frame, 'Enter csv filename for XY table', pth, '',
        'XY table file (*.csv)|*.csv',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
    try:
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            filename = os.path.splitext(filename)[0]+'.csv'
            File = open(filename,'w')
        else:
            filename = None
    finally:
        dlg.Destroy()
    if not filename:
        return
    for i in range(len(XY)):
        if len(names):
            header = '%s,%s(%s)\n'%(labelX,labelY,names[i])
        else:
            header = '%s,%s(%d)\n'%(labelX,labelY,i)
        File.write(header)
        for x,y in XY[i].T:
            File.write('%.3f,%.3f\n'%(x,y))
    File.close()
    print (' XY data saved to: '+filename)

def PeakListSave(G2frame,file,peaks):
    'Save powder peaks to a data file'
    print ('save peak list to file: '+G2frame.peaklistfile)
    if not peaks:
        dlg = wx.MessageDialog(G2frame, 'No peaks!', 'Nothing to save!', wx.OK)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()
        return
    for peak in peaks:
        file.write("%10.4f %12.2f %10.3f %10.3f \n" % \
            (peak[0],peak[2],peak[4],peak[6]))
    print ('peak list saved')

def IndexPeakListSave(G2frame,peaks):
    'Save powder peaks from the indexing list'
    file = open(G2frame.peaklistfile,'wa')
    print ('save index peak list to file: '+G2frame.peaklistfile)
    wx.BeginBusyCursor()
    try:
        if not peaks:
            dlg = wx.MessageDialog(G2frame, 'No peaks!', 'Nothing to save!', wx.OK)
            try:
                dlg.ShowModal()
            finally:
                dlg.Destroy()
            return
        for peak in peaks:
            file.write("%12.6f\n" % (peak[7]))
        file.close()
    finally:
        wx.EndBusyCursor()
    print ('index peak list saved')

######################################################################
def ExportPowderList(G2frame):
    '''Returns a list of extensions supported by :func:`ExportPowder`
    along with their descriptions (note that a extension may be repeated
    but descriptions are unique).
    This is used in :meth:`GSASIIimgGUI.AutoIntFrame` only.

    :param wx.Frame G2frame: the GSAS-II main data tree window
    '''
    extList = []
    extLabel = []
    for obj in G2frame.exporterlist:
        if 'powder' in obj.exporttype:
            try:
                obj.Writer
                extList.append(obj.extension)
                extLabel.append(obj.formatName)
            except AttributeError:
                pass
    return extList,extLabel

def ExportPowder(G2frame,TreeName,fileroot,extension,hint=''):
    '''Writes a single powder histogram using the Export routines.
    This is used in :meth:`GSASIIimgGUI.AutoIntFrame` only.

    :param wx.Frame G2frame: the GSAS-II main data tree window
    :param str TreeName: the name of the histogram (PWDR ...) in the data tree
    :param str fileroot: name for file to be written, extension ignored
    :param str extension: extension for file to be written (start with '.'). Must
      match a powder export routine that has a Writer object.
    :param str hint: a string that must match the export's format
    '''
    filename = os.path.abspath(os.path.splitext(fileroot)[0]+extension)
    for obj in G2frame.exporterlist:
        if obj.extension == extension and 'powder' in obj.exporttype:
            if hint and hint not in obj.formatName: continue
            obj.currentExportType = 'powder'
            obj.InitExport(None)
            obj.loadTree() # load all histograms in tree into dicts
            if TreeName not in obj.Histograms:
                raise Exception('Histogram not found: '+str(TreeName))
            try:
                obj.Writer
            except AttributeError:
                continue
            try:
                obj.Writer(TreeName,filename)
                print('wrote file '+filename)
                return
            except Exception:
                print('Export Routine for '+extension+' failed.')
    else:
        print('No Export routine supports extension '+extension)

def ExportSequentialFullCIF(G2frame,seqData,Controls):
    '''Handles access to CIF exporter a bit differently for sequential fits, as this is
    not accessed via the usual export menus
    '''
    from exports import G2export_CIF
    ##################### debug code to reload exporter before each use ####
    #import importlib as imp
    #imp.reload(G2export_CIF)
    #print('reload G2export_CIF')
    ########################################################################
    obj = G2export_CIF.ExportProjectCIF(G2frame)
    obj.Exporter(None,seqData=seqData,Controls=Controls)

def ExportSequential(G2frame,data,obj,exporttype):
    '''
    Used to export from every phase/dataset in a sequential refinement using
    a .Writer method for either projects or phases. Prompts to select histograms
    and for phase exports, which phase(s).

    :param wx.Frame G2frame: the GSAS-II main data tree window
    :param dict data: the sequential refinement data object
    :param exporter obj: an exporter object
    :param str exporttype: indicates the type of export ('project' or 'phase')
    '''
    if len(data['histNames']) == 0:
        G2G.G2MessageBox(G2frame,'There are no sequential histograms','Warning')
    obj.InitExport(None)
    obj.loadTree()
    obj.loadParmDict()
    if len(data['histNames']) == 1:
        histlist = data['histNames']
    else:
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select histograms to export from list',
                                 'Select histograms',data['histNames'])
        if dlg.ShowModal() == wx.ID_OK:
            histlist = [data['histNames'][l] for l in dlg.GetSelections()]
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
    if exporttype == 'Phase':
        phaselist = list(obj.Phases.keys())
        if len(obj.Phases) == 0:
            G2G.G2MessageBox(G2frame,'There are no phases in sequential ref.','Warning')
            return
        elif len(obj.Phases) > 1:
            dlg = G2G.G2MultiChoiceDialog(G2frame,'Select phases to export from list',
                                    'Select phases', phaselist)
            if dlg.ShowModal() == wx.ID_OK:
                phaselist = [phaselist[l] for l in dlg.GetSelections()]
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
        filename = obj.askSaveFile()
        if not filename: return True
        obj.dirname,obj.filename = os.path.split(filename)
        print('Writing output to file '+str(obj.filename)+"...")
        mode = 'w'
        for p in phaselist:
            for h in histlist:
                obj.SetSeqRef(data,h)
                #GSASIIpath.IPyBreak()
                obj.Writer(h,phasenam=p,mode=mode)
                mode = 'a'
        print('...done')
    elif exporttype == 'Project':  # note that the CIF exporter is not yet ready for this
        filename = obj.askSaveFile()
        if not filename: return True
        obj.dirname,obj.filename = os.path.split(filename)
        print('Writing output to file '+str(obj.filename)+"...")
        mode = 'w'
        for h in histlist:
            obj.SetSeqRef(data,h)
            obj.Writer(h,mode=mode)
            print('\t'+str(h)+' written')
            mode = 'a'
        print('...done')
    elif exporttype == 'Powder':
        filename = obj.askSaveFile()
        if not filename: return True
        obj.dirname,obj.filename = os.path.split(filename)
        print('Writing output to file '+str(obj.filename)+"...")
        mode = 'w'
        for h in histlist:
            obj.SetSeqRef(data,h)
            obj.Writer(h,mode=mode)
            print('\t'+str(h)+' written')
            mode = 'a'
        print('...done')

def ReadDIFFaX(DIFFaXfile):
    print ('read '+DIFFaXfile)
    Layer = {'Laue':'-1','Cell':[False,1.,1.,1.,90.,90.,90,1.],'Width':[[10.,10.],[False,False]],
        'Layers':[],'Stacking':[],'Transitions':[],'Toler':0.01,'AtInfo':{}}
    df = open(DIFFaXfile,'r')
    lines = df.readlines()
    df.close()
    struct = False
    Struct = []
    stack = False
    Stack = []
    trans = False
    Trans = []
    for diff in lines:
        diff = diff[:-1].lower()
        if '!'  in diff:
            continue
        while '}' in diff: #strip comments
            iB = diff.index('{')
            iF = diff.index('}')+1
            if iB:
                diff = diff[:iB]
            else:
                diff = diff[iF:]
        if not diff:
            continue
        if diff.strip() == 'instrumental':
            continue
        if diff.strip() == 'structural':
            struct = True
            continue
        elif diff.strip() == 'stacking':
            struct = False
            stack = True
            continue
        elif diff.strip() == 'transitions':
            stack = False
            trans = True
            continue
        diff = diff.strip()
        if struct:
            if diff:
                Struct.append(diff)
        elif stack:
            if diff:
                Stack.append(diff)
        elif trans:
            if diff:
                Trans.append(diff)

#STRUCTURE records
    laueRec = Struct[1].split()
    Layer['Laue'] = laueRec[0]
    if Layer['Laue'] == 'unknown' and len(laueRec) > 1:
        Layer['Toler'] = float(laueRec[1])    #tolerance for 'unknown'?
    if Layer['Laue'] == '2/m(1)': Layer['Laue'] = '2/m(c)'
    if Layer['Laue'] == '2/m(2)': Layer['Laue'] = '2/m(ab)'
    cell = Struct[0].split()
    Layer['Cell'] = [False,float(cell[0]),float(cell[1]),float(cell[2]),90.,90.,float(cell[3]),1.0]
    nLayers = int(Struct[2])
    N = 3
    if 'layer' not in Struct[3]:
        N = 4
        if Struct[3] != 'infinite':
            width = Struct[3].split()
            Layer['Width'][0] = [float(width[0]),float(width[1])]
    for nL in range(nLayers):
        if '=' in Struct[N]:
            name = Struct[N].split('=')
            sameas = int(name[1])-1
            Layer['Layers'].append({'Name':name[0],'SameAs':Layer['Layers'][sameas]['Name'],'Symm':'None','Atoms':[]})
            N += 1
            continue
        Symm = 'None'
        if 'centro' in Struct[N+1]: Symm = '-1'
        Layer['Layers'].append({'Name':Struct[N],'SameAs':'','Symm':Symm,'Atoms':[]})
        N += 2
        while 'layer' not in Struct[N]:
            atom = Struct[N][4:].split()
            atomType = G2el.FixValence(Struct[N][:4].replace(' ','').strip().capitalize())
            if atomType not in Layer['AtInfo']:
                Layer['AtInfo'][atomType] = G2el.GetAtomInfo(atomType)
            atomName = '%s(%s)'%(atomType,atom[0])
            newVals = []
            for val in atom[1:6]:
                if '/' in val:
                    newVals.append(eval(val+'.'))
                else:
                    newVals.append(float(val))
            atomRec = [atomName,atomType,newVals[0],newVals[1],newVals[2],newVals[4],newVals[3]/78.9568]
            Layer['Layers'][-1]['Atoms'].append(atomRec)
            N += 1
            if N > len(Struct)-1:
                break
#TRANSITIONS records
    transArray = []
    N = 0
    for i in range(nLayers):
        transArray.append([])
        for j in range(nLayers):
            vals = Trans[N].split()
            newVals = []
            for val in vals[:4]:
                if '/' in val:
                    newVals.append(eval(val+'.'))
                else:
                    newVals.append(float(val))
            transArray[-1].append(newVals+['',False])
            N += 1
    Layer['Transitions'] = transArray
#STACKING records
    Layer['Stacking'] = [Stack[0],'']
    if Stack[0] == 'recursive':
        Layer['Stacking'][1] = Stack[1]
    elif Stack[0] == 'explicit':
        if Stack[1] == 'random':
            Layer['Stacking'][1] = Stack[1]
        else:
            Layer['Stacking'][1] = 'list'
            Layer['Stacking'].append('')
            for stack in Stack[2:]:
                Layer['Stacking'][2] += ' '+stack
    return Layer

def saveNewPhase(G2frame,phData,newData,phlbl,msgs,orgFilName):
    '''create a .gpx file from a structure from the BilbaoSite pseudosym site
    saved in newData
    '''
    def fmtCell(cell):
        s = ''
        for i in cell[0:3]: s += f"{i:.3f}, "
        for i in cell[3:5]: s += f"{i:.2f}, "
        s += f"{cell[5]:.2f}"
        return s
    if newData is None:
        print(phlbl,'empty structure')
        return
    elif type(newData) is str:
        msgs[phlbl] = newData
        return
    # create a new phase
    try:
        sgnum = int(newData[0].strip())
        sgsym = G2spc.spgbyNum[sgnum]
        sgname = sgsym.replace(" ","")
    except:
        print(f'Problem with processing record:\n{newData}')
        return
    newPhase = copy.deepcopy(phData)
    newPhase['ranId'] = ran.randint(0,sys.maxsize),
    if 'magPhases' in phData: del newPhase['magPhases']
    generalData = newPhase['General']
    generalData['SGData'] = SGData = G2spc.SpcGroup(sgsym)[1]
    generalData['Cell'][1:7] = [float(i) for i in newData[1].split()]
    generalData['Cell'][7] = G2lat.calc_V(G2lat.cell2A(generalData['Cell'][1:7]))
    cx,ct,cs,cia = generalData['AtomPtrs']
    Atoms = newPhase['Atoms'] = []
    for a in newData[3:]:
        if not a.strip(): continue
        try:
            elem,n,wyc,x,y,z = a.split()
            atom = []
            atom.append(elem+n)
            atom.append(elem)
            atom.append('')
            for i in x,y,z: atom.append(float(i))
            atom.append(1.0)
            SytSym,Mult = G2spc.SytSym(np.array(atom[3:6]),SGData)[:2]
            atom.append(SytSym)
            atom.append(Mult)
            atom.append('I')
            atom += [0.02,0.,0.,0.,0.,0.,0.,]
            atom.append(ran.randint(0,sys.maxsize))
            Atoms.append(atom)
        except:
            print(f'error in atom line {a}')
        #finally: pass
    phData.update(newPhase)
    G2elem.SetupGeneral(phData,phData['General']['Mydir'])  # fixup composition info
    # save new file
    G2frame.GSASprojectfile = os.path.splitext(orgFilName
                            )[0]+'_super_'+sgname.replace('/','$')+'.gpx'
    while os.path.exists(G2frame.GSASprojectfile):
        s = re.split(r'_([\d]+)\.gpx',G2frame.GSASprojectfile)
        if len(s) == 1:
            G2frame.GSASprojectfile = os.path.splitext(G2frame.GSASprojectfile)[0] + '_1.gpx'
        else:
            num = 10
            try:
                num = int(s[1]) + 1
            except:
                pass
            G2frame.GSASprojectfile = f'{s[0]}_{num}.gpx'
    ProjFileSave(G2frame)
    # get transformed contents
    nacomp,nccomp = G2mth.phaseContents(phData)
    msgs[phlbl] = f"With space group {sgsym} and cell={fmtCell(generalData['Cell'][1:7])}"
    msgs[phlbl] += f", vol={generalData['Cell'][7]:.2f} A^3"
    msgs[phlbl] += f", project file created as {G2frame.GSASprojectfile}"

    msgs[phlbl] += f". After transform, unit cell {G2mth.fmtPhaseContents(nccomp)}"
    msgs[phlbl] += f", density={G2mth.getDensity(generalData)[0]:.2f} g/cm^3"
    msgs[phlbl] += f". Asymmetric unit {G2mth.fmtPhaseContents(nacomp)} ({len(phData['Atoms'])} atoms)."
    return G2frame.GSASprojectfile

def mkParmDictfromTree(G2frame,sigDict=None):
    '''Load the GSAS-II refinable parameters from the tree into dict parmDict
    Updating refined values to those from the last cycle. Optionally
    compute the s.u. values for the parameters and place them in sigDict.

    The actions in the routine are used in a number of places around
    the GSAS-II code where it would be "cleaner" to use this instead.
    Perhaps that will happen as code revisions are made. One example
    of this, :meth:`GSASIIfiles.ExportBaseclass.loadParmDict`.

    :param wx.Frame G2frame: a reference to the main GSAS-II window
    :param dict sigDict: a Python dict with sigma (s.u.) values for
      each parameter

    :returns: parmDict, a dict with the value for all refined and most
      unrefined GSAS-II parameters used in the diffraction computations.
      This parmDict version has only values, as opposed to the version
      used in some parts of the code that has refinement flags and initial
      values as well.
    '''
    from . import GSASIIstrIO as G2stIO
    from . import GSASIIstrMath as G2stMth
    from . import GSASIImapvars as G2mv
    G2frame.CheckNotebook()
    parmDict = {}
    rigidbodyDict = {}
    covDict = {}
    consDict = {}

    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    if G2frame.GPXtree.IsEmpty(): return # nothing to do
    rigidbodyDict = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Rigid bodies'))
    covDict = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Covariance'))
    consDict = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Constraints'))

    rbVary,rbDict =  G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
    parmDict.update(rbDict)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,maxSSwave =  \
        G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False)
    parmDict.update(phaseDict)
    hapVary,hapDict,controlDict =  G2stIO.GetHistogramPhaseData(Phases,Histograms,Print=False,resetRefList=False)
    parmDict.update(hapDict)
    histVary,histDict,controlDict =  G2stIO.GetHistogramData(Histograms,Print=False)
    parmDict.update(histDict)
    parmDict.update(zip(
            covDict.get('varyList',[]),
            covDict.get('variables',[])))
    if sigDict is not None:
        sigDict.update(dict(zip(
            covDict.get('varyList',[]),
            covDict.get('sig',[])
            )))

    # expand to include constraints: first compile a list of constraints
    constList = []
    for item in consDict:
        if item.startswith('_'): continue
        constList += consDict[item]
    G2mv.InitVars()     # process constraints
    constrDict,fixedList,ignored = G2mv.ProcessConstraints(constList)
    varyList = list(covDict.get('varyListStart',[]))
    if varyList is None and len(constrDict) == 0:
        # no constraints can use varyList
        varyList = covDict.get('varyList')
    elif varyList is None:
        varyList = []
        # # old GPX file from before pre-constraint varyList is saved
        print (' *** Old refinement: Please use Calculate/Refine to redo  ***')
        return None
    varyList = list(varyList)
    # add constraint generated parameters to parmDict
    G2mv.EvaluateMultipliers(constrDict,parmDict)
    errmsg,warnmsg,groups,parmlist = G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict)
    # enforce constraints on values (should have been done; no effect)
    G2mv.Map2Dict(parmDict,varyList)   # N.B. changes varyList
    G2mv.Dict2Map(parmDict)
    # add uncertainties for constrained & RB parameters
    if sigDict is not None:
        sigDict.update(G2mv.ComputeDepESD(
            covDict['covMatrix'],covDict['varyList']))
        sigDict.update(G2stMth.computeRBsu(parmDict,Phases,rigidbodyDict,
            covDict['covMatrix'],covDict['varyList'],covDict['sig']))
    return parmDict

def LogCellChanges(G2frame):
    '''Log varied cell parameters into the data tree notebook'''
    Controls = G2frame.GPXtree.GetItemPyData(
        G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
    if Controls['max cyc'] == 0: return # no fit so no change
    covData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Covariance'))
    parmDict = mkParmDictfromTree(G2frame)
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    for phase in Phases:
        phasedict = Phases[phase]
        pId = phasedict['pId']
        SGData = phasedict['General']['SGData']
        for hist in phasedict['Histograms']:
            if not phasedict['Histograms'][hist]['Use']: continue
            #skip single crystal histograms!
            if 'Type' in phasedict['Histograms'][hist] and 'S' in phasedict['Histograms'][hist]['Type']: continue
            hId = Histograms[hist]['hId']
            if any(phasedict['Histograms'][hist]['HStrain'][1]) or phasedict['General']['Cell'][0]:
                cellList,cellSig = G2lat.getCellSU(pId,hId,SGData,parmDict,covData)
                txt = G2lat.showCellSU(cellList,cellSig,SGData)
                G2frame.AddToNotebook(f'Phase {pId} Hist {hId}: {txt}',
                                          'CEL',False)
if __name__ == '__main__':
    from . import GSASIIdataGUI
    application = GSASIIdataGUI.GSASIImain(0)
    G2frame = application.main
    #app = wx.PySimpleApp()
    #G2frame = wx.Frame(None) # create a frame
    #frm.Show(True)
    #filename = '/tmp/notzip.zip'
    #filename = '/tmp/all.zip'
    #filename = '/tmp/11bmb_7652.zip'

    #selection=None, confirmoverwrite=True, parent=None
    # choicelist=[ ('a','b','c'),
    #              ('test1','test2'),('no choice',)]
    # titles = [ 'a, b or c', 'tests', 'No option here']
    # dlg = MultipleChoicesDialog(
    #     choicelist,titles,
    #     parent=frm)
    # if dlg.ShowModal() == wx.ID_OK:
    #     print 'Got OK'
    #imagefile = '/tmp/NDC5_00237_3.ge3'
    #Comments, Data, Npix, Image = G2fil.GetImageData(G2frame,imagefile,imageOnly=False,ImageTag=None)

    #print("\n\nResults loaded to Comments, Data, Npix and Image\n\n")

    #GSASIIpath.IPyBreak_base()
