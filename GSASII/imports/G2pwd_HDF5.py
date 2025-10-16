# -*- coding: utf-8 -*-
'''Use to read powder patterns from HDF5 files. At present the only supported 
format is a NeXus variant named NXazint1d. 
'''

from __future__ import division, print_function
import os

try:
    import h5py
except ImportError:
    h5py = None
import numpy as np
from .. import GSASIIobj as G2obj
from .. import GSASIIfiles as G2fil

#from .. import GSASIIpath
#breakpoint = GSASIIpath.IPyBreak_base

class HDF5_Reader(G2obj.ImportPowderData):
    '''Routine to read multiple powder patterns from an HDF5 file. 

    This importer targets NXazint1d and NXazint2d NeXus files from 
    MAX IV. 
    Perhaps in the future, other types of HDF5 powder data sources as well. 

    The main file is <file>.hdf or <file>.h5, but optionally sample and 
    instrument parameters can be placed in <file>.samprm and <file>.instprm. 
    Any parameters placed in that file will override values set in the HDF5
    file. 
    '''
    #mode = None
    def __init__(self):
        if h5py is None:
            self.UseReader = False
            msg = 'HDF5 Reader skipped because h5py module is not installed'
            G2fil.ImportErrorMsg(msg,{'HDF5 importer':['h5py','hdf5']})
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hdf','.h5'),strictExtension=True,
            formatName = 'MAX IV HDF5',longFormatName = 'MaxIV NXazint1d HDF5 integrated scans')
        self.scriptable = True
        #self.Iparm = {} #only filled for EDS data

    def ShowH5Element(self,obj,keylist):
        '''Format the contents of an HDF5 entry as a single line. Not used for 
        reading files, only used in :meth:`HDF5list` which is here for software
        development. 
        '''
        k = '/'.join(keylist)
        l = obj.get(k, getlink=True)
        if isinstance(l, h5py.ExternalLink): 
            return f'link to file {l.filename}'
        try:
            typ = str(type(obj[k]))
        except:
            return f'**Error** with key {k}'
            
        if ".Dataset'" in typ:
            datfmt = obj[k].dtype
            if datfmt == 'O' or str(datfmt).startswith('|S'):
                # byte string
                return f'value={obj[k][()].decode()}'
            elif datfmt == 'bool': # Bool
                return f'value={bool(obj[k][()])}'
            elif datfmt in ('<f8', 'uint8', 'int64', '<f4'): # scalar value or array of values
                try:
                    len(obj[k][()])
                    return f'array {obj[k].shape}'
                except:
                    return f'value={obj[k][()]}'
            else:
                return f'dataset of type {repr(datfmt)}'
        elif ".Group'" in typ:
            return "(group)"
        else:
            return f'type is {type(obj[k])}'

    def RecurseH5Element(self,obj,prefix=[],length=None):
        '''Returns a list of entries of all keys in the HDF5 file
        (or group) in `obj`. Note that `obj` can be a file object, created by 
        `h5py.File` or can be a subset `fp['key/subkey']`.
        
        If length is specified, only the entries with arrays of that
        length are returned.

        The returned list is organized where: 
          * entry 0 is the top-level keys (/a, /b,...),
          * entry 1 has the first level keys (/a/c /a/d, /b/d, /b/e,...)
          * ...
        '''
        try:
            self.HDF5entries
        except AttributeError:
            self.HDF5entries = []
        depth = len(prefix)
        if len(self.HDF5entries) < depth+1:
            self.HDF5entries.append([])
        for i in obj:
            nextprefix = prefix+[i]
            if length is None:
                self.HDF5entries[depth].append(nextprefix)
            try:
                typ = str(type(obj[i]))
            except:
                print(f'**Error** with key {prefix}/{i}')
                continue
            if length is not None and ".Group'" not in typ:
                # get length of this obj[i]
                try:
                    if len(obj[i]) == length:
                        self.HDF5entries[depth].append(nextprefix)
                except TypeError:
                    continue
            # check for link objects
            l = obj.get(i, getlink=True)
            if isinstance(l, h5py.ExternalLink): continue
            if ".Group'" in typ:
                #t = f'{prefix}/{i}'
                #print(f'\n{nextprefix} contents {(60-len(t))*'='}')
                self.RecurseH5Element(obj[i],nextprefix,length)
        return self.HDF5entries
        
                
    def HDF5list(self, filename):
        '''Shows the contents of an HDF5 file as a short listing. 
        This is not used for HDF5 reading, but is of help with a new
        type of HDF5 file to see what is present.

        :param filename: 
        '''
        def ShowH5NeXusName(obj,keylist):
            key = '/'.join(keylist)
            if "NX_class" in obj[key].attrs:
                return obj[key].attrs["NX_class"]

        fp = h5py.File(filename, 'r')
        #print(f'Contents of {filename}')
        HDF5entries = self.RecurseH5Element(fp)
        strings = []
        for i,j in enumerate(HDF5entries):
            if not strings or strings[-1] != 60*'=': 
                strings.append(60*'=')
            m = 0
            for k in j:
                m = max(m,len('/'.join(k)))
            for k in j:
                nxname = ShowH5NeXusName(fp,k)
                lbl = self.ShowH5Element(fp,k)
                if '\n' in lbl:
                    lbl = '; '.join(lbl.split('\n'))
                if len(lbl) > 50:
                    lbl = lbl[:50] + '...'
                # if '\n' in lbl:
                #     lbl = lbl.split()[0] + '...'
                if lbl != '(group)': strings.append(f"{'/'.join(k):{m}s} {lbl}")
                if nxname: print(f"{'/'.join(k):{m}s} {lbl} {nxname}")
        with open(filename+'_contents.txt', 'w') as fp:
            for i in strings: fp.write(f'{i}\n')
                    
    def ContentsValidator(self, filename):
        '''Test if valid by seeing if the HDF5 library recognizes the file. 
        Then get file type (currently MAX IV NeXus/NXazint[12]d only)
        '''
        try:
            fp = h5py.File(filename, 'r')
            # test for MaxIV NeXus/NXazint1d & NXazint2d
            test = True
            while test:
                test = False
                entry = getNeXusBase(fp)
                if entry is None: break # not NeXus
                if 'definition' not in fp[entry]: break # not MaxIV NXazint*
                definition = fp[entry+'/definition'][()].decode()
                if definition == 'NXazint1d': return True
                if definition == 'NXazint2d': return True
            # test for next HDF5 type here
            #
        except IOError: # not HDF5
            return False
        finally:
            fp.close()
        return False        # nothing passed -- not valid

    def Reader(self, filename, ParentFrame=None, **kwarg):
        '''Scan file for sections needed by defined file types (currently 
        MAX IV NeXus/NXazint[12]d only) 
        and then use appropriate routine to read the file.

        Since usually there will be lots of scans in a single file, 
        the goal is that the first pass should read the file into 
        a buffer (if available) and subsequent calls will not 
        need to access the file. 
        '''
        fpbuffer = kwarg.get('buffer',{})
        if not hasattr(self,'blknum'):
            if self.selections is None or len(self.selections) == 0:
                self.blknum = 0
            else:
                self.blknum = min(self.selections)
        # was file read into buffer? If so skip opening file to save time
        definition = fpbuffer.get('definition','')
        if definition == 'NXazint1d':
            return self.readNXazint1d(filename, fpbuffer)
        elif definition == 'NXazint2d':
            return self.readNXazint2d(filename, fpbuffer)
        # first or non-buffered read
        try:
            fp = h5py.File(filename, 'r')
            entry = getNeXusBase(fp)
            if entry: # NeXus
                if 'definition' in fp[entry]: # MaxIV NXazint*
                    definition = fp[entry+'/definition'][()].decode()
                    fpbuffer['definition'] = definition
                    if definition == 'NXazint1d':
                        return self.readNXazint1d(filename, fpbuffer)
                    elif definition == 'NXazint2d':
                        return self.readNXazint2d(filename, fpbuffer)
            # not a supported file type
            return False
        except IOError:
            print ('cannot open file '+ filename)
            return False
        finally:
            fp.close()

        print (f'Unknown type of HDF5 powder file {filename}')
        return False

    def readNXazint1d(self, filename, fpbuffer={}):
        '''Read HDF5 file in NeXus as produced by MAX IV as a NXazint1d
        see https://nxazint-hdf5-nexus-3229ecbd09ba8a773fbbd8beb72cace6216dfd5063e1.gitlab-pages.esrf.fr/classes/contributed_definitions/NXazint1d.html
        '''
        #self.instmsg = 'HDF file'
        self.comments = []
        doread = False # has the file already been read into a buffer?
        fileItems = {
            # arrays
            'radial_axis':('NXdata','radial_axis'),
            'I':('NXdata','I'),
            'I_errors':('NXdata','I_errors'),
            # floats
            'wavelength':('NXmonochromator','wavelength'),
            'polarization_factor':('NXparameters','polarization_factor'),
            # strings
            'instrument/name':('NXinstrument','name'),
            'unit':('NXparameters','unit'),
            'sample/name':('NXsample','name'),
            'source/name':('NXsource','name'),
        }
        # test if we have what we need in the buffer
        for k in fileItems:
            if k not in fpbuffer:                
                doread = True
                break
        if doread:
            # Nope, need to fill the buffer
            try:
                fp = h5py.File(filename, 'r')
                entry = getNeXusBase(fp)
                # lookup NeXus locations
                nexusDict = {i:None for i in set([i[0] for i in fileItems.values()])}
                recurseNeXusEntries(fp,entry,nexusDict)
                # save selected items from file in buffer
                savedKeys = []
                for k,loc in fileItems.items():
                    if nexusDict[loc[0]] is None:
                        fpbuffer[k] = None
                        continue
                    key = '/'.join((nexusDict[loc[0]],)+loc[1:])
                    savedKeys.append(key)
                    if key not in fp:
                        fpbuffer[k] = None
                        continue
                    val = fp[key]
                    if val.shape:
                        fpbuffer[k] = np.array(val)
                    elif 'float' in str(val.dtype):
                        fpbuffer[k] = float(val[()])
                        self.comments.append(f'{k}={val[()]}')
                    else:
                        fpbuffer[k] = val[()].decode()
                        self.comments.append(f'{k}={fpbuffer[k]}')
                self.numbanks = len(fpbuffer['I'])
                # save arrays that are potentially tracking the parametric conditions
                # into ParamTrackingVars.
                # e.g. variables with the same length as the humber of datasets
                fpbuffer['ParamTrackingVars'] = {}
                paramItems = self.RecurseH5Element(fp,length=self.numbanks)
                for i in paramItems:
                    for j in i:
                        key = '/'.join(j)
                        if key in savedKeys: continue # standard data array
                        obj = fp.get(key)
                        if obj is None: continue
                        if len(obj[()].shape) != 1: continue
                        # are all values the same? If so, put them into the comments
                        # for the first histogram only. If they are changing, note that 
                        # here and later they will be put into every histogram.
                        if all(obj[0] == obj):
                            self.comments.append(f'{key.split("/")[-1]}={obj[0]}')
                        else:
                            fpbuffer['ParamTrackingVars'][key] = np.array(obj[()])
            except IOError:
                print (f'Can not open or read file {filename}')
                return False
            finally:
                fp.close()
            if fpbuffer['unit'] != '2th':
                print('NXazint1d HDF5 file has units',fpbuffer['entry/reduction/input/unit'])
                self.errors = 'NXazint1d only can be read with 2th units'
                return False
            # initialize the block selection
            if self.selections is None or len(self.selections) == 0:
                self.blknum = 0
            else:
                self.blknum = min(self.selections)
        # now pull the selected dataset from the buffer
        x = fpbuffer['radial_axis']
        y = fpbuffer['I'][self.blknum]
        try:
            esd = fpbuffer['I_errors'][self.blknum]
            w = np.where(esd==0,0,np.nan_to_num(1/esd**2))
        except:
            w = np.nan_to_num(1/y)    # best we can do, alas w/o reported s.u.'s
        self.powderdata = [x,y,w,np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)]
        # add parametric var as a comment
        for key,arr in fpbuffer['ParamTrackingVars'].items():
            val = arr[self.blknum]
            self.comments.append(f'{key.split("/")[-1]}={val}')
            if 'temperature' in key:
                self.Sample['Temperature'] = val # in K already
            elif 'time' in key:
                self.Sample['Time'] = val # should be seconds
            elif 'chi' in key:
                self.Sample['Chi'] = val # not sure if correct mapping
            elif 'phi' in key:
                self.Sample['Phi'] = val
            elif 'omega' in key:
                self.Sample['Omega'] = val
        self.powderentry[0] = filename
        #self.powderentry[1] = Pos # position offset (never used, I hope)
        self.powderentry[2] = self.blknum  # bank number
        self.idstring = f'#{self.blknum} {os.path.split(filename)[1][:60]}'
        self.instdict['wave'] = fpbuffer['wavelength']
        # if not, are there more [selected] images that after this to be read?
        self.repeat = False
        if self.blknum < self.numbanks-1:
            if self.selections is None or len(self.selections) == 0:
                self.blknum += 1
                self.repeat = True
            else:
                try:
                    s = sorted(self.selections)
                    self.blknum = s[s.index(self.blknum)+1]
                    self.repeat = True
                except IndexError:   # last selected image has been read
                    self.repeat = False
        return True

# NeXus support routines. These were influenced heavily by Frederik Holm Gjørup
# Also see NeXus support in plaid (https://github.com/fgjorup/plaid/blob/main/plaid/nexus.py)

def getNeXusBase(fp):
    '''This returns the base entry in a NeXus compilant HDF5 file
    (usually "/entry" for MaxIV files) or None if this is not a valid 
    NeXus file. 
    '''
    for key in fp:
        if ("NX_class" in fp[key].attrs and
                fp[key].attrs["NX_class"] == "NXentry"):
            return key

def getNeXusEntry(fp,base,target):
    '''This returns the entry in a NeXus compilant HDF5 file matching 
    the name target, or None, if this is not found as a child of the key `base`.
    Not in use as it is more practical to use :func:`recurseNeXusEntries`.
    '''
    for key in fp[base]:
        subkey = '/'.join([base,key])
        if "NX_class" in fp[subkey].attrs:
            #print(key, list(fp[subkey].attrs),fp[subkey].attrs["NX_class"])
            if ("NX_class" in fp[subkey].attrs and
                    fp[subkey].attrs["NX_class"] == target):
                return subkey
        else:
            print(key)

def recurseNeXusEntry(fp,node,target):
    '''Recurse through the HDF5 tree looking for NeXus class `target`. 
    Not in use, as :func:`recurseNeXusEntries` is used to get all 
    targets in a single pass through the tree.
    '''
    if node is None: return  # needed?
    val = fp[node]
    if ("NX_class" in val.attrs and
                val.attrs["NX_class"] == target):
        return node
    if not isinstance(val, h5py.Group): return
    for key in val:
        subkey = '/'.join([node,key])
        res = recurseNeXusEntry(fp,subkey,target)
        if res: return res

def recurseNeXusEntries(fp,node,targetdict):
    '''recurse through the HDF5 tree looking for the NeXus classes
    in `targetdict`, storing the HDF5 key for each class in the dict

    :param fp: HDF5 file pointer
    :param str node: name of current HDF5 key
    :param dict targetdict: dict to place HDF5 keys corresponding to
       the desired NeXus classes. As input this has the NeXus classes
       is the dict keys and the 
    '''
    val = fp[node]
    if ("NX_class" in val.attrs and val.attrs["NX_class"] in targetdict):
        targetdict[val.attrs["NX_class"]] = node
    if isinstance(val, h5py.Group): 
        for key in val:
            recurseNeXusEntries(fp,'/'.join([node,key]),targetdict)
