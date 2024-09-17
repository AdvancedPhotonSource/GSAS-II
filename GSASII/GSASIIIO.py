# -*- coding: utf-8 -*-
'''
Misc routines for input and output, including image reading follow. 

TODO: This module needs some work to separate wx from non-wx routines. GUI 
routines should probably move to GSASIIctrlGUI.
'''

from __future__ import division, print_function

# Allow this to be imported without wx present. Was needed for G2scriptable, but is
# likely not needed anymore
try:
    import wx
except ImportError:
    # was needed by sphinx, but probably not anymore
    class Placeholder(object):
        def __init__(self):
            self.Dialog = object
    wx = Placeholder()
import math
import numpy as np
import numpy.ma as ma

import copy
import platform
import pickle as cPickle
import sys
import re
import random as ran
import GSASIIpath
try:
    import GSASIIdataGUI as G2gd
except ImportError:
    pass
import GSASIIobj as G2obj
import GSASIIlattice as G2lat
try:
    import GSASIIpwdGUI as G2pdG
    import GSASIIimgGUI as G2imG
except ImportError:
    pass
import GSASIIElem as G2el
import GSASIIstrIO as G2stIO
import GSASIImapvars as G2mv
import GSASIIfiles as G2fil
try:
    import GSASIIctrlGUI as G2G
except ImportError:
    pass
import os
import os.path as ospath

DEBUG = False       #=True for various prints
TRANSP = False      #=true to transpose images for testing
if GSASIIpath.GetConfigValue('Transpose'): TRANSP = True
npsind = lambda x: np.sin(x*np.pi/180.)

def sfloat(S):
    'Convert a string to float. An empty field or a unconvertable value is treated as zero'
    if S.strip():
        try:
            return float(S)
        except ValueError:
            pass
    return 0.0

def sint(S):
    'Convert a string to int. An empty field is treated as zero'
    if S.strip():
        return int(S)
    else:
        return 0

def trim(val):
    '''Simplify a string containing leading and trailing spaces
    as well as newlines, tabs, repeated spaces etc. into a shorter and
    more simple string, by replacing all ranges of whitespace
    characters with a single space. 

    :param str val: the string to be simplified

    :returns: the (usually) shortened version of the string
    '''
    return re.sub('\\s+', ' ', val).strip()

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
            dlg.SetFilename(''+ospath.split(imagefile)[1])
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
    if ImageTag:
        TreeLbl += ' #'+'%04d'%(ImageTag)
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

def GetImageData(G2frame,imagefile,imageOnly=False,ImageTag=None,FormatName=''):
    '''Read a single image with an image importer. This is called to reread an image
    after it has already been imported with :meth:`GSASIIdataGUI.GSASII.OnImportGeneric`
    (or :func:`ReadImages` in Auto Integration) so it is not necessary to reload metadata.

    :param wx.Frame G2frame: main GSAS-II Frame and data object.
    :param str imagefile: name of image file
    :param bool imageOnly: If True return only the image,
      otherwise  (default) return more (see below)
    :param int/str ImageTag: specifies a particular image to be read from a file.
      First image is read if None (default).
    :param str formatName: the image reader formatName

    :returns: an image as a numpy array or a list of four items:
      Comments, Data, Npix and the Image, as selected by imageOnly

    '''
    # determine which formats are compatible with this file
    primaryReaders = []
    secondaryReaders = []
    for rd in G2frame.ImportImageReaderlist:
        flag = rd.ExtensionValidator(imagefile)
        if flag is None: 
            secondaryReaders.append(rd)
        elif flag:
            if not FormatName:
                primaryReaders.append(rd)
            elif FormatName in rd.formatName:       #This is a kluge because the rd.formatName was changed!
                primaryReaders.append(rd)
    if len(secondaryReaders) + len(primaryReaders) == 0:
        print('Error: No matching format for file '+imagefile)
        raise Exception('No image read')
    errorReport = ''
    if not imagefile:
        return
    for rd in primaryReaders+secondaryReaders:
        rd.ReInitialize() # purge anything from a previous read
        rd.errors = "" # clear out any old errors
        if not rd.ContentsValidator(imagefile): # rejected on cursory check
            errorReport += "\n  "+rd.formatName + ' validator error'
            if rd.errors: 
                errorReport += ': '+rd.errors
                continue
        if imageOnly:
            ParentFrame = None # prevent GUI access on reread
        else:
            ParentFrame = G2frame
        if GSASIIpath.GetConfigValue('debug'):
            flag = rd.Reader(imagefile,ParentFrame,blocknum=ImageTag)
        else:
            flag = False
            try:
                flag = rd.Reader(imagefile,ParentFrame,blocknum=ImageTag)
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
            #rd.readfilename = imagefile
            if imageOnly:
                return rd.Image
            else:
                return rd.Comments,rd.Data,rd.Npix,rd.Image
    else:
        print('Error reading file '+imagefile)
        print('Error messages(s)\n'+errorReport)
        raise Exception('No image read')    

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
    cPickle.dump([Comments,Data,Npix,image],File,2)
    File.close()
    return

objectScanIgnore = [int,bool,float,str,np.float64,np.float32,np.int32,np.int64,np.int16,np.ndarray,G2obj.G2VarObj,G2obj.ExpressionObj,np.bool_]
try:
    objectScanIgnore += [ma.MaskedArray] # fails in doc builds
except AttributeError:
    pass
    
if '2' in platform.python_version_tuple()[0]:
    objectScanIgnore += [unicode,long,]
    
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
    
def cPickleLoad(fp):
    if '2' in platform.python_version_tuple()[0]:
        return cPickle.load(fp)
    else:
       return cPickle.load(fp,encoding='latin-1')

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
                data = cPickleLoad(fp) # first block in file should be Phases
                if data[0][0] != 'Phases':
                    raise Exception('Unexpected block in {} file. How did this happen?'
                            .format(GPXphase))
                Phases = {}
                for name,vals in data[1:]:
                    Phases[name] = vals
                name,CovData = cPickleLoad(fp)[0] # 2nd block in file should be Covariance
                name,RigidBodies = cPickleLoad(fp)[0] # 3rd block in file should be Rigid Bodies
                fp.close()
                # index the histogram updates
                hist = open(GPXhist,'rb')
                try:
                    while True:
                        loc = hist.tell()
                        datum = cPickleLoad(hist)[0]
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
                data = cPickleLoad(filep)
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
                hdata = cPickleLoad(hist)
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
            G2G.updateNotifier(G2frame,int(LastSavedUsing))
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
        # (Controls should always be created at this point)
        try:
            Controls = G2frame.GPXtree.GetItemPyData(
                G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
            Controls['LastSavedAs'] = os.path.abspath(G2frame.GSASprojectfile)
            Controls['LastSavedUsing'] = str(GSASIIpath.GetVersionNumber())
            Controls['PythonVersions'] = G2frame.PackageVersions
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
                cPickle.dump(data,file,2)
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
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Substances'),G2pdG.SetDefaultSubstances())
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Sample Parameters'),Sample)
            G2frame.GPXtree.SetItemPyData(G2frame.GPXtree.AppendItem(Id,text='Models'),G2pdG.SetDefaultSASDModel())
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

def striphist(var,insChar=''):
    'strip a histogram number from a var name'
    sv = var.split(':')
    if len(sv) <= 1: return var
    if sv[1]:
        sv[1] = insChar
    return ':'.join(sv)

######################################################################
# base classes for reading various types of data files
#   not used directly, only by subclassing
######################################################################
class ExportBaseclass(object):
    '''Defines a base class for the exporting of GSAS-II results.

    This class is subclassed in the various exports/G2export_*.py files. Those files
    are imported in :meth:`GSASIIdataGUI.GSASII._init_Exports` which defines the
    appropriate menu items for each one and the .Exporter method is called
    directly from the menu item.

    Routines may also define a .Writer method, which is used to write a single
    file without invoking any GUI objects.
    '''
    # TODO: review exporters producing exceptions where .Writer can't be used where G2frame is None (see CIF)
    # TODO: review conflicting uses of .Writer with mode (SeqRef) & elsewhere
    # TODO: move this class to G2fil
    def __init__(self,G2frame,formatName,extension,longFormatName=None,):
        self.G2frame = G2frame
        self.formatName = formatName # short string naming file type
        self.extension = extension
        if longFormatName: # longer string naming file type
            self.longFormatName = longFormatName
        else:
            self.longFormatName = formatName
        self.OverallParms = {}
        self.Phases = {}
        self.Histograms = {}
        self.powderDict = {}
        self.sasdDict = {}
        self.refdDict = {}
        self.xtalDict = {}
        self.parmDict = {}
        self.sigDict = {}
        self.fp = None
        # updated in InitExport:
        self.currentExportType = None # type of export that has been requested
        # updated in ExportSelect (when used):
        self.phasenam = None # a list of selected phases
        self.histnam = None # a list of selected histograms
        self.filename = None # name of file to be written (single export) or template (multiple files)
        self.dirname = '' # name of directory where file(s) will be written
        self.fullpath = '' # name of file being written -- full path
        
        # items that should be defined in a subclass of this class
        self.exporttype = []  # defines the type(s) of exports that the class can handle.
        # The following types are defined: 'project', "phase", "powder", "single"
        self.multiple = False # set as True if the class can export multiple phases or histograms
        # self.multiple is ignored for "project" exports

    def InitExport(self,event):
        '''Determines the type of menu that called the Exporter and
        misc initialization. 
        '''
        self.filename = None # name of file to be written (single export)
        self.dirname = '' # name of file to be written (multiple export)
        if event:
            self.currentExportType = self.G2frame.ExportLookup.get(event.Id)

    def MakePWDRfilename(self,hist):
        '''Make a filename root (no extension) from a PWDR histogram name

        :param str hist: the histogram name in data tree (starts with "PWDR ")
        '''
        file0 = ''
        file1 = hist[5:]
        # replace repeated blanks
        while file1 != file0:
            file0 = file1
            file1 = file0.replace('  ',' ').strip()
        file0 = file1.replace('Azm= ','A')
        # if angle has unneeded decimal places on aziumuth, remove them
        if file0[-3:] == '.00': file0 = file0[:-3]
        file0 = file0.replace('.','_')
        file0 = file0.replace(' ','_')
        return file0

    def ExportSelect(self,AskFile='ask'):
        '''Selects histograms or phases when needed. Sets a default file name when
        requested into self.filename; always sets a default directory in self.dirname.

        :param bool AskFile: Determines how this routine processes getting a
          location to store the current export(s).
          
          * if AskFile is 'ask' (default option), get the name of the file to be written;
            self.filename and self.dirname are always set. In the case where
            multiple files must be generated, the export routine should do this
            based on self.filename as a template.
          * if AskFile is 'dir', get the name of the directory to be used;
            self.filename is not used, but self.dirname is always set. The export routine
            will always generate the file name.
          * if AskFile is 'single', get only the name of the directory to be used when
            multiple items will be written (as multiple files) are used
            *or* a complete file name is requested when a single file
            name is selected. self.dirname is always set and self.filename used
            only when a single file is selected.  
          * if AskFile is 'default', creates a name of the file to be used from
            the name of the project (.gpx) file. If the project has not been saved,
            then the name of file is requested.
            self.filename and self.dirname are always set. In the case where
            multiple file names must be generated, the export routine should do this
            based on self.filename.
          * if AskFile is 'default-dir', sets self.dirname from the project (.gpx)
            file. If the project has not been saved, then a directory is requested.
            self.filename is not used.

        :returns: True in case of an error
        '''
        
        numselected = 1
        if self.currentExportType == 'phase':
            if len(self.Phases) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any phases.')
                return True
            elif len(self.Phases) == 1:
                self.phasenam = list(self.Phases.keys())
            elif self.multiple: 
                choices = sorted(self.Phases.keys())
                phasenum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if phasenum is None: return True
                self.phasenam = [choices[i] for i in phasenum]
                if not self.phasenam: return True
                numselected = len(self.phasenam)
            else:
                choices = sorted(self.Phases.keys())
                phasenum = G2G.ItemSelector(choices,self.G2frame)
                if phasenum is None: return True
                self.phasenam = [choices[phasenum]]
                numselected = len(self.phasenam)
        elif self.currentExportType == 'single':
            if len(self.xtalDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any single crystal data.')
                return True
            elif len(self.xtalDict) == 1:
                self.histnam = list(self.xtalDict.values())
            elif self.multiple:
                choices = sorted(self.xtalDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.xtalDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'powder':
            if len(self.powderDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any powder data.')
                return True
            elif len(self.powderDict) == 1:
                self.histnam = list(self.powderDict.values())
            elif self.multiple:
                choices = sorted(self.powderDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.powderDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'sasd':
            if len(self.sasdDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any small angle data.')
                return True
            elif len(self.sasdDict) == 1:
                self.histnam = self.sasdDict.values()
            elif self.multiple:
                choices = sorted(self.sasdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.sasdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'refd':
            if len(self.refdDict) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any reflectivity data.')
                return True
            elif len(self.refdDict) == 1:
                self.histnam = self.refdDict.values()
            elif self.multiple:
                choices = sorted(self.refdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=True)
                if not hnum: return True
                self.histnam = [choices[i] for i in hnum]
                numselected = len(self.histnam)
            else:
                choices = sorted(self.refdDict.values())
                hnum = G2G.ItemSelector(choices,self.G2frame)
                if hnum is None: return True
                self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        elif self.currentExportType == 'image':
            if len(self.Histograms) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any images.')
                return True
            elif len(self.Histograms) == 1:
                self.histnam = list(self.Histograms.keys())
            else:
                choices = sorted(self.Histograms.keys())
                hnum = G2G.ItemSelector(choices,self.G2frame,multiple=self.multiple)
                if self.multiple:
                    if not hnum: return True
                    self.histnam = [choices[i] for i in hnum]
                else:
                    if hnum is None: return True
                    self.histnam = [choices[hnum]]
                numselected = len(self.histnam)
        if self.currentExportType == 'map':
            # search for phases with maps
            mapPhases = []
            choices = []
            for phasenam in sorted(self.Phases):
                phasedict = self.Phases[phasenam] # pointer to current phase info            
                if len(phasedict['General']['Map'].get('rho',[])):
                    mapPhases.append(phasenam)
                    if phasedict['General']['Map'].get('Flip'):
                        choices.append('Charge flip map: '+str(phasenam))
                    elif phasedict['General']['Map'].get('MapType'):
                        choices.append(
                            str(phasedict['General']['Map'].get('MapType'))
                            + ' map: ' + str(phasenam))
                    else:
                        choices.append('unknown map: '+str(phasenam))
            # select a map if needed
            if len(mapPhases) == 0:
                self.G2frame.ErrorDialog(
                    'Empty project',
                    'Project does not contain any maps.')
                return True
            elif len(mapPhases) == 1:
                self.phasenam = mapPhases
            else: 
                phasenum = G2G.ItemSelector(choices,self.G2frame,multiple=self.multiple)
                if self.multiple:
                    if not phasenum: return True
                    self.phasenam = [mapPhases[i] for i in phasenum]
                else:
                    if phasenum is None: return True
                    self.phasenam = [mapPhases[phasenum]]
            numselected = len(self.phasenam)

        # items selected, now set self.dirname and usually self.filename 
        if AskFile == 'ask' or (AskFile == 'single' and numselected == 1) or (
            AskFile == 'default' and not self.G2frame.GSASprojectfile
            ):
            filename = self.askSaveFile()
            if not filename: return True
            self.dirname,self.filename = os.path.split(filename)
        elif AskFile == 'dir' or AskFile == 'single' or (
            AskFile == 'default-dir' and not self.G2frame.GSASprojectfile
            ):
            self.dirname = self.askSaveDirectory()
            if not self.dirname: return True
        elif AskFile == 'default-dir' or AskFile == 'default':
            self.dirname,self.filename = os.path.split(
                os.path.splitext(self.G2frame.GSASprojectfile)[0] + self.extension
                )
        else:
            raise Exception('This should not happen!')

    def loadParmDict(self,computeSU=False):
        '''Load the GSAS-II refinable parameters from the tree into a dict (self.parmDict). Update
        refined values to those from the last cycle and set the uncertainties for the
        refined parameters in another dict (self.sigDict).

        Expands the parm & sig dicts to include values derived from constraints.

        This could be made faster for sequential fits as info for each histogram is loaded
        later when iterating over histograms. 
        '''
        self.G2frame.CheckNotebook()
        self.parmDict = {}
        self.sigDict = {} # dict with s.u. values, currently used only for CIF & Bracket exports
        self.RBsuDict = {} # dict with s.u. values for atoms in a rigid body
        self.constList = [] # constraint entries from data tree
        rigidbodyDict = {}
        covDict = {}
        consDict = {}
        Histograms,Phases = self.G2frame.GetUsedHistogramsAndPhasesfromTree()
        if self.G2frame.GPXtree.IsEmpty(): return # nothing to do
        item, cookie = self.G2frame.GPXtree.GetFirstChild(self.G2frame.root)
        while item:
            name = self.G2frame.GPXtree.GetItemText(item)
            if name == 'Rigid bodies':
                 rigidbodyDict = self.G2frame.GPXtree.GetItemPyData(item)
            elif name == 'Covariance':
                 covDict = self.G2frame.GPXtree.GetItemPyData(item)
            elif name == 'Constraints':
                 consDict = self.G2frame.GPXtree.GetItemPyData(item)
            item, cookie = self.G2frame.GPXtree.GetNextChild(self.G2frame.root, cookie)
        rbVary,rbDict =  G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)
        self.parmDict.update(rbDict)
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,maxSSwave =  \
            G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False)
        self.parmDict.update(phaseDict)
        hapVary,hapDict,controlDict =  G2stIO.GetHistogramPhaseData(Phases,Histograms,Print=False,resetRefList=False)
        self.parmDict.update(hapDict)
        histVary,histDict,controlDict =  G2stIO.GetHistogramData(Histograms,Print=False)
        self.parmDict.update(histDict)
        self.parmDict.update(zip(
            covDict.get('varyList',[]),
            covDict.get('variables',[])))
        self.sigDict = dict(zip(
            covDict.get('varyList',[]),
            covDict.get('sig',[])))
        # expand to include constraints: first compile a list of constraints
        self.constList = []
        for item in consDict:
            if item.startswith('_'): continue
            self.constList += consDict[item]
        # now process the constraints
        G2mv.InitVars()
        constrDict,fixedList,ignored = G2mv.ProcessConstraints(self.constList)
        varyList = covDict.get('varyListStart')
        if varyList is None and len(constrDict) == 0:
            # no constraints can use varyList
            varyList = covDict.get('varyList')
        elif varyList is None:
            varyList = []
            # # old GPX file from before pre-constraint varyList is saved
            # print (' *** Old refinement: Please use Calculate/Refine to redo  ***')
            # raise Exception(' *** Export aborted ***')
        else:
            varyList = list(varyList)
        # add symmetry-generated constraints
        rigidbodyDict = self.G2frame.GPXtree.GetItemPyData(   
            G2gd.GetGPXtreeItemId(self.G2frame,self.G2frame.root,'Rigid bodies'))
        rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
        rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,Print=False)  # done twice, needed?
        Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,EFtables,ORBtables,BLtables,MFtables,maxSSwave = \
            G2stIO.GetPhaseData(Phases,RestraintDict=None,rbIds=rbIds,Print=False) # generates atom symmetry constraints
        msg = G2mv.EvaluateMultipliers(constrDict,phaseDict)
        if msg:
            print('Unable to interpret multiplier(s): '+msg)
            raise Exception(' *** CIF creation aborted ***')
        errmsg,warnmsg,groups,parmlist = G2mv.GenerateConstraints(varyList,constrDict,fixedList,self.parmDict)
        if errmsg:
            # this really should not happen
            print (' *** ERROR - constraints are internally inconsistent ***')
            print ('Errors: ',errmsg)
            if warnmsg: print ('Warnings'+warnmsg)
            raise Exception(' *** CIF creation aborted ***')
        G2mv.Map2Dict(self.parmDict,varyList)   # changes varyList
        G2mv.Dict2Map(self.parmDict)   # add the constrained values to the parameter dictionary
        # and add their uncertainties into the esd dictionary (sigDict)
        if covDict.get('covMatrix') is not None and computeSU:
            self.sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],covDict['varyList'],noSym=True))
            if 'depSigDict' in self.OverallParms['Covariance']:
                self.sigDict.update(
                    {i:v[1] for i,v in self.OverallParms['Covariance']['depSigDict'].items()})
            # compute the s.u.'s on rigid bodies
            import GSASIIstrMath as G2stMth
            self.RBsuDict = G2stMth.computeRBsu(self.parmDict,Phases,rigidbodyDict,
                            covDict['covMatrix'],covDict['varyList'],covDict['sig'])

    def loadTree(self):
        '''Load the contents of the data tree into a set of dicts
        (self.OverallParms, self.Phases and self.Histogram as well as self.powderDict
        & self.xtalDict)
        
        * The childrenless data tree items are overall parameters/controls for the
          entire project and are placed in self.OverallParms
        * Phase items are placed in self.Phases
        * Data items are placed in self.Histogram. The key for these data items
          begin with a keyword, such as PWDR, IMG, HKLF,... that identifies the data type.
        '''
        self.OverallParms = {}
        self.powderDict = {}
        self.sasdDict = {}
        self.refdDict = {}
        self.xtalDict = {}
        self.Phases = {}
        self.Histograms = {}
        self.SeqRefdata = None
        self.SeqRefhist = None
        self.DelayOpen = False
        if self.G2frame.GPXtree.IsEmpty(): return # nothing to do
        histType = None        
        if self.currentExportType == 'phase':
            # if exporting phases load them here
            sub = G2gd.GetGPXtreeItemId(self.G2frame,self.G2frame.root,'Phases')
            if not sub:
                print ('no phases found')
                return True
            item, cookie = self.G2frame.GPXtree.GetFirstChild(sub)
            while item:
                phaseName = self.G2frame.GPXtree.GetItemText(item)
                self.Phases[phaseName] =  self.G2frame.GPXtree.GetItemPyData(item)
                item, cookie = self.G2frame.GPXtree.GetNextChild(sub, cookie)
            # Get rigid body info into self.OverallParms 
            for key in ('Rigid bodies','Covariance'):
                item = G2gd.GetGPXtreeItemId(self.G2frame,self.G2frame.root,key)
                if item:
                    self.OverallParms[key] = self.G2frame.GPXtree.GetItemPyData(item)
                item, cookie = self.G2frame.GPXtree.GetNextChild(sub, cookie)
            return
        elif self.currentExportType == 'single':
            histType = 'HKLF'
        elif self.currentExportType == 'powder':
            histType = 'PWDR'
        elif self.currentExportType == 'image':
            histType = 'IMG'
        elif self.currentExportType == 'sasd':
            histType = 'SASD'
        elif self.currentExportType == 'refd':
            histType = 'REFD'

        if histType: # Loading just one kind of tree entry
            item, cookie = self.G2frame.GPXtree.GetFirstChild(self.G2frame.root)
            while item:
                name = self.G2frame.GPXtree.GetItemText(item)
                if name.startswith(histType):
                    if self.Histograms.get(name): # there is already an item with this name
                        print('Histogram name '+str(name)+' is repeated. Renaming')
                        if name[-1] == '9':
                            name = name[:-1] + '10'
                        elif name[-1] in '012345678':
                            name = name[:-1] + str(int(name[-1])+1)
                        else:                            
                            name += '-1'
                    self.Histograms[name] = {}
                    # the main info goes into Data, but the 0th
                    # element contains refinement results, carry
                    # that over too now. 
                    self.Histograms[name]['Data'] = self.G2frame.GPXtree.GetItemPyData(item)[1]
                    self.Histograms[name][0] = self.G2frame.GPXtree.GetItemPyData(item)[0]
                    item2, cookie2 = self.G2frame.GPXtree.GetFirstChild(item)
                    while item2: 
                        child = self.G2frame.GPXtree.GetItemText(item2)
                        self.Histograms[name][child] = self.G2frame.GPXtree.GetItemPyData(item2)
                        item2, cookie2 = self.G2frame.GPXtree.GetNextChild(item, cookie2)
                item, cookie = self.G2frame.GPXtree.GetNextChild(self.G2frame.root, cookie)
            # index powder and single crystal histograms by number
            for hist in self.Histograms:
                if hist.startswith("PWDR"): 
                    d = self.powderDict
                elif hist.startswith("HKLF"): 
                    d = self.xtalDict
                elif hist.startswith("SASD"):
                    d = self.sasdDict
                elif hist.startswith("REFD"):
                    d = self.refdDict
                else:
                    return                    
                i = self.Histograms[hist].get('hId')
                if i is None and not d.keys():
                    i = 0
                elif i is None or i in d.keys():
                    i = max(d.keys())+1
                d[i] = hist
            return
        # else standard load: using all interlinked phases and histograms
        self.Histograms,self.Phases = self.G2frame.GetUsedHistogramsAndPhasesfromTree()
        item, cookie = self.G2frame.GPXtree.GetFirstChild(self.G2frame.root)
        while item:
            name = self.G2frame.GPXtree.GetItemText(item)
            if name == 'Restraints':
                self.OverallParms[name] = self.G2frame.GPXtree.GetItemPyData(item)
            else:
                item2, cookie2 = self.G2frame.GPXtree.GetFirstChild(item)
                if not item2: 
                    self.OverallParms[name] = self.G2frame.GPXtree.GetItemPyData(item)
            item, cookie = self.G2frame.GPXtree.GetNextChild(self.G2frame.root, cookie)
        # index powder and single crystal histograms
        for hist in self.Histograms:
            i = self.Histograms[hist]['hId']
            if hist.startswith("PWDR"): 
                self.powderDict[i] = hist
            elif hist.startswith("HKLF"): 
                self.xtalDict[i] = hist
            elif hist.startswith("SASD"):
                self.sasdDict[i] = hist
            elif hist.startswith("REFD"):
                self.refdDict[i] = hist

    def dumpTree(self,mode='type'):
        '''Print out information on the data tree dicts loaded in loadTree.
        Used for testing only.
        '''
        if self.SeqRefdata and self.SeqRefhist:
            print('Note that dumpTree does not show sequential results')
        print ('\nOverall')
        if mode == 'type':
            def Show(arg): return type(arg)
        else:
            def Show(arg): return arg
        for key in self.OverallParms:
            print ('  '+key+Show(self.OverallParms[key]))
        print ('Phases')
        for key1 in self.Phases:
            print ('    '+key1+Show(self.Phases[key1]))
        print ('Histogram')
        for key1 in self.Histograms:
            print ('    '+key1+Show(self.Histograms[key1]))
            for key2 in self.Histograms[key1]:
                print ('      '+key2+Show(self.Histograms[key1][key2]))

    def defaultSaveFile(self):
        return os.path.abspath(
            os.path.splitext(self.G2frame.GSASprojectfile
                             )[0]+self.extension)
        
    def askSaveFile(self):
        '''Ask the user to supply a file name

        :returns: a file name (str) or None if Cancel is pressed

        '''
        #pth = G2G.GetExportPath(self.G2frame)
        if self.G2frame.GSASprojectfile:
            defnam = os.path.splitext(
                os.path.split(self.G2frame.GSASprojectfile)[1]
                )[0]+self.extension
        else:
            defnam = 'default' + self.extension
        return G2G.askSaveFile(self.G2frame,defnam,self.extension,self.longFormatName)

    def askSaveDirectory(self):
        '''Ask the user to supply a directory name. Path name is used as the
        starting point for the next export path search. 

        :returns: a directory name (str) or None if Cancel is pressed

        TODO: Can this be replaced with G2G.askSaveDirectory?
        '''
        pth = G2G.GetExportPath(self.G2frame)
        dlg = wx.DirDialog(
            self.G2frame, 'Input directory where file(s) will be written', pth,
            wx.DD_DEFAULT_STYLE)
        dlg.CenterOnParent()
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                self.G2frame.LastExportDir = filename
            else:
                filename = None
        finally:
            dlg.Destroy()
        return filename

    # Tools for file writing. 
    def OpenFile(self,fil=None,mode='w',delayOpen=False):
        '''Open the output file

        :param str fil: The name of the file to open. If None (default)
          the name defaults to self.dirname + self.filename.
          If an extension is supplied, it is not overridded,
          but if not, the default extension is used. 
        :param str mode:  The mode can 'w' to write a file, or 'a' to append to it. If 
          the mode is 'd' (for debug), output is displayed on the console. 
        :returns: the file object opened by the routine which is also
          saved as self.fp
        '''
        if mode == 'd': # debug mode
            self.fullpath = '(stdout)'
            self.fp = sys.stdout
            return
        if not fil:
            if not os.path.splitext(self.filename)[1]:
                self.filename += self.extension
            fil = os.path.join(self.dirname,self.filename)
        self.fullpath = os.path.abspath(fil)
        self.DelayOpen = False
        if delayOpen:
            self.DelayOpen = True
            self.fp = None
            return
        self.fp = open(self.fullpath,mode)
        return self.fp

    def openDelayed(self,mode='w'):
        self.DelayOpen = False
        self.fp = open(self.fullpath,mode)
        return self.fp

    def Write(self,line):
        '''write a line of output, attaching a line-end character

        :param str line: the text to be written. 
        '''
        if self.fp is None:
            raise Exception('Attempt to Write without use of OpenFile')
        self.fp.write(line+'\n')
        
    def CloseFile(self,fp=None):
        '''Close a file opened in OpenFile

        :param file fp: the file object to be closed. If None (default)
          file object self.fp is closed. 
        '''
        if self.fp is None and self.DelayOpen:
            if GSASIIpath.GetConfigValue('debug'): 
                print('Delayed open: close before file not created')
            return
        if self.fp is None:
            if GSASIIpath.GetConfigValue('debug'): 
                raise Exception('Attempt to CloseFile without use of OpenFile')
            else:
                print('Attempt to CloseFile without use of OpenFile')
                return
        if self.fp == sys.stdout: return # debug mode
        if fp is None:
            fp = self.fp
            self.fp = None
        if fp is not None: fp.close()
        
    def SetSeqRef(self,data,hist):
        '''Set the exporter to retrieve results from a sequential refinement
        rather than the main tree
        '''
        self.SeqRefdata = data
        self.SeqRefhist = hist
        data_name = data[hist]
        for i,val in zip(data_name['varyList'],data_name['sig']):
            self.sigDict[i] = val
            self.sigDict[striphist(i)] = val
        for i in data_name['parmDict']:
            self.parmDict[striphist(i)] = data_name['parmDict'][i]
            self.parmDict[i] = data_name['parmDict'][i]
            # zero out the dA[xyz] terms, they would only bring confusion
            key = i.split(':')
            if len(key) < 3: continue
            if key[2].startswith('dA'):
                self.parmDict[i] = 0.0
        for i,(val,sig) in data_name.get('depParmDict',{}).items():
            self.parmDict[i] = val
            self.sigDict[i] = sig
        #GSASIIpath.IPyBreak()

    # Tools to pull information out of the data arrays
    def GetCell(self,phasenam,unique=False):
        """Gets the unit cell parameters and their s.u.'s for a selected phase

        :param str phasenam: the name for the selected phase
        :param bool unique: when True, only directly refined parameters
          (a in cubic, a & alpha in rhombohedral cells) are assigned 
          positive s.u. values. Used as True for CIF generation.
        :returns: `cellList,cellSig` where each is a 7 element list corresponding
          to a, b, c, alpha, beta, gamma, volume where `cellList` has the
          cell values and `cellSig` has their uncertainties.
        """
        if self.SeqRefdata and self.SeqRefhist:
            return self.GetSeqCell(phasenam,self.SeqRefdata[self.SeqRefhist])
        phasedict = self.Phases[phasenam] # pointer to current phase info
        try:
            pfx = str(phasedict['pId'])+'::'
            A,sigA = G2stIO.cellFill(pfx,phasedict['General']['SGData'],self.parmDict,self.sigDict)
            cellSig = G2stIO.getCellEsd(pfx,phasedict['General']['SGData'],A,
                self.OverallParms['Covariance'],unique=unique)  # returns 7 vals, includes sigVol
            cellList = G2lat.A2cell(A) + (G2lat.calc_V(A),)
            return cellList,cellSig
        except KeyError:
            cell = phasedict['General']['Cell'][1:]
            return cell,7*[0]
            
    def GetSeqCell(self,phasenam,data_name):
        """Gets the unit cell parameters and their s.u.'s for a selected phase
        and histogram in a sequential fit

        :param str phasenam: the name for the selected phase
        :param dict data_name: the sequential refinement parameters for the selected histogram
        :returns: `cellList,cellSig` where each is a 7 element list corresponding
          to a, b, c, alpha, beta, gamma, volume where `cellList` has the
          cell values and `cellSig` has their uncertainties.
        """
        phasedict = self.Phases[phasenam]
        SGdata = phasedict['General']['SGData']
        pId = phasedict['pId']
        RecpCellTerms = G2lat.cell2A(phasedict['General']['Cell'][1:7])
        ESDlookup = {}
        Dlookup = {}
        varied = [striphist(i) for i in data_name['varyList']]
        for item,val in data_name['newCellDict'].items():
            if item in varied:
                ESDlookup[val[0]] = item
                Dlookup[item] = val[0]
        A = RecpCellTerms[:]
        for i in range(6):
            var = str(pId)+'::A'+str(i)
            if var in ESDlookup:
                A[i] = data_name['newCellDict'][ESDlookup[var]][1] # override with refined value
        cellDict = dict(zip([str(pId)+'::A'+str(i) for i in range(6)],A))
        zeroDict = {i:0.0 for i in cellDict}
        A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata,cellDict,zeroDict)
        covData = {
            'varyList': [Dlookup.get(striphist(v),v) for v in data_name['varyList']],
            'covMatrix': data_name['covMatrix']
            }
        return list(G2lat.A2cell(A)) + [G2lat.calc_V(A)], G2stIO.getCellEsd(str(pId)+'::',SGdata,A,covData)
                
    def GetAtoms(self,phasenam):
        """Gets the atoms associated with a phase. Can be used with standard
        or macromolecular phases

        :param str phasenam: the name for the selected phase
        :returns: a list of items for eac atom where each item is a list containing:
          label, typ, mult, xyz, and td, where

          * label and typ are the atom label and the scattering factor type (str)
          * mult is the site multiplicity (int)
          * xyz is contains a list with four pairs of numbers:
            x, y, z and fractional occupancy and
            their standard uncertainty (or a negative value)
          * td is contains a list with either one or six pairs of numbers:
            if one number it is U\\ :sub:`iso` and with six numbers it is
            U\\ :sub:`11`, U\\ :sub:`22`, U\\ :sub:`33`, U\\ :sub:`12`, U\\ :sub:`13` & U\\ :sub:`23`
            paired with their standard uncertainty (or a negative value)
        """
        phasedict = self.Phases[phasenam] # pointer to current phase info            
        cx,ct,cs,cia = phasedict['General']['AtomPtrs']
        cfrac = cx+3
        fpfx = str(phasedict['pId'])+'::Afrac:'        
        atomslist = []
        for i,at in enumerate(phasedict['Atoms']):
            if phasedict['General']['Type'] == 'macromolecular':
                label = '%s_%s_%s_%s'%(at[ct-1],at[ct-3],at[ct-4],at[ct-2])
            else:
                label = at[ct-1]
            fval = self.parmDict.get(fpfx+str(i),at[cfrac])
            fsig = self.sigDict.get(fpfx+str(i),-0.009)
            mult = at[cs+1]
            typ = at[ct]
            xyz = []
            for j,v in enumerate(('x','y','z')):
                val = at[cx+j]
                pfx = str(phasedict['pId']) + '::A' + v + ':' + str(i)
                val = self.parmDict.get(pfx, val)
                dpfx = str(phasedict['pId'])+'::dA'+v+':'+str(i)
                sig = self.sigDict.get(dpfx,-0.000009)
                xyz.append((val,sig))
            xyz.append((fval,fsig))
            td = []
            if at[cia] == 'I':
                pfx = str(phasedict['pId'])+'::AUiso:'+str(i)
                val = self.parmDict.get(pfx,at[cia+1])
                sig = self.sigDict.get(pfx,-0.0009)
                td.append((val,sig))
            else:
                for i,var in enumerate(('AU11','AU22','AU33','AU12','AU13','AU23')):
                    pfx = str(phasedict['pId'])+'::'+var+':'+str(i)
                    val = self.parmDict.get(pfx,at[cia+2+i])
                    sig = self.sigDict.get(pfx,-0.0009)
                    td.append((val,sig))
            atomslist.append((label,typ,mult,xyz,td))
        return atomslist
######################################################################
def ExportPowderList(G2frame):
    '''Returns a list of extensions supported by :func:`GSASIIIO:ExportPowder`
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
    import G2export_CIF
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

def postURL(URL,postdict,getcookie=None,usecookie=None,
                timeout=None,retry=2,mode='get'):
    '''Posts a set of values as from a web form using the "get" or "post" 
    protocols. 
    If access fails to an https site, the access is retried with http.

    :param str URL: the URL to post; typically something 
       like 'http://www.../dir/page?'
    :param dict postdict: contains keywords and values, such
       as {'centrosymmetry': '0', 'crystalsystem': '0', ...}
    :param dict getcookie: dict to save cookies created in call, or None
       (default) if not needed. 
    :param dict usecookie: dict containing cookies to be used in call, 
       or None (default) if not needed.
    :param int timeout: specifies a timeout period for the get or post (default 
      is None, which means the timeout period is set by the server). The value 
      when specified is the time in seconds to wait before giving up on the 
      request.
    :param int retry: the number of times to retry the request, if it times out.
      This is only used if timeout is specified. The default is 2. Note that
      if retry is left at the default value (2), The timeout is increased by
      25% for the second try.
    :param str mode: either 'get' (default) or 'post'. Determines how
       the request will be submitted. 
    :returns: a string with the response from the web server or None
       if access fails.
    '''
    try:
        import requests # delay this until now, since rarely needed
    except:
        # this import seems to fail with the Anaconda pythonw on
        # macs; it should not!
        print('Warning: failed to import requests. Python config error')
        return None

    if mode == 'get':
        reqopt = requests.get
    else:
        reqopt = requests.post
    
    repeat = True
    count = 0
    while repeat:
        count += 1
        r = None
        repeat = False
        try:
            if timeout is not None:
                r = reqopt(URL,params=postdict,cookies=usecookie,
                                     timeout=timeout)
            else:
                r = reqopt(URL,params=postdict,cookies=usecookie)
            if r.status_code == 200:
                if GSASIIpath.GetConfigValue('debug'): print('request OK')
                page = r.text
                if getcookie is not None:
                    getcookie.update(r.cookies)
                return page # success
            else:
                print('request to {} failed. Reason={}'.format(URL,r.reason))
        except requests.exceptions.ConnectionError as msg:
            if 'time' in str(msg) and 'out' in str(msg): 
                print(f'server timeout accessing {URL}')
                if GSASIIpath.GetConfigValue('debug'): print('full error=',msg)
                if timeout is not None and count < retry:
                    if retry == 2:
                        timeout *= 1.25
                        print(f'retry with timout={timeout} sec')
                    repeat = True
            else:
                print('connection error - not on internet?')
                if URL.startswith('https:'):
                    print('Retry with http://')
                    repeat = True
                    URL = URL.replace('https:','http:')
        except requests.exceptions.Timeout as msg:
            print(f'timeout accessing {URL}')
            if GSASIIpath.GetConfigValue('debug'): print('full error=',msg)
            if timeout is not None and count < retry:
                if retry == 2:
                    timeout *= 1.25
                    print(f'retry with timout={timeout} sec')
                repeat = True
            if timeout is not None and count <= retry: repeat = True
        except requests.exceptions.ReadTimeout as msg:
            print(f'timeout reading from {URL}')
            if GSASIIpath.GetConfigValue('debug'): print('full error=',msg)
            if timeout is not None and count < retry:
                if retry == 2:
                    timeout *= 1.25
                    print(f'retry with timout={timeout} sec')
                repeat = True
        except requests.exceptions.ConnectTimeout:
            print(f'timeout accessing {URL}')
        except Exception as msg:    # other error
            print(f'Error accessing {URL}')
            if GSASIIpath.GetConfigValue('debug'): print(msg)
        finally:
            if r: r.close()
    else:
        return None

if __name__ == '__main__':
    import GSASIIdataGUI
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
    imagefile = '/tmp/NDC5_00237_3.ge3'
    Comments, Data, Npix, Image = GetImageData(G2frame,imagefile,imageOnly=False,ImageTag=None)

    print("\n\nResults loaded to Comments, Data, Npix and Image\n\n")

    GSASIIpath.IPyBreak_base()
