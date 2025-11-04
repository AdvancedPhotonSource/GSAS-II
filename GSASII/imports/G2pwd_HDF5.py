# -*- coding: utf-8 -*-
'''Use to read powder patterns from HDF5 files. At present the only supported 
format are two NeXus variants from MaxIV named NXazint1d and NXazint2d,
but this can be expanded to handle more HDF5/NeXus formats
'''
import os
try:
    import h5py
except ImportError:
    h5py = None
import numpy as np
from .. import GSASIIobj as G2obj
from .. import GSASIIfiles as G2fil

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
    def __init__(self):
        if h5py is None:
            self.UseReader = False
            msg = 'HDF5 Reader skipped because h5py module is not installed'
            G2fil.ImportErrorMsg(msg,{'HDF5 importer':['h5py','hdf5']})
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.hdf','.h5'),strictExtension=True,
            formatName = 'MAXIV NeXus',longFormatName = 'Max IV NXazintXd NeXus integrated scans')
        self.scriptable = True

    def ContentsValidator(self, filename):
        '''Test if valid by seeing if the HDF5 library recognizes the file. 
        Then get file type (currently MAX IV NXazint[12]d (NeXus) only)
        '''
        try:
            definition = ''
            fp = h5py.File(filename, 'r')
            # test for MaxIV NeXus/NXazint1d & NXazint2d
            test = True
            while test: # block for standard NXazint1d and NXazint2d files,
                        # use break to bail out and try next block
                test = False
                entry = getNeXusBase(fp)
                if entry is None: break # not NeXus
                if 'definition' not in fp[entry]: break # not MaxIV NXazint*
                definition = fp[entry+'/definition'][()].decode()
                # get names for datasets so we can select them
                if definition == 'NXazint1d':
                    fileItems = {
                        'I':('NXdata','I'),
                        'unit':('NXparameters','unit'),
                        }
                    buffer = {}
                    if not self.readInNeXus(filename,buffer,fileItems,'NXazint1d-validate'):
                        return False
                    nhist = len(buffer['I'])
                    self.selections = list(range(nhist))
                    for i in range(nhist):
                        self.dnames.append(f'#{i} {os.path.split(filename)[1][:60]}')
                    return True
                if definition == 'NXazint2d':
                    fileItems = {
                        'I':('NXdata','I'),
                        'unit':('NXparameters','unit'),
                        'azimuthal_axis':('NXdata','azimuthal_axis'),
                    }
                    buffer = {}
                    if not self.readInNeXus(filename,buffer,fileItems,'NXazint2d-validate'):
                        return False
                    #numazimuth = buffer['azimuth_bins']
                    numazimuth = len(buffer['azimuthal_axis'])
                    numbanks = len(buffer['I'])
                    nhist = numbanks * numazimuth
                    self.selections = list(range(nhist))
                    for i in range(nhist):
                        # group by parametric variable
                        numScan = i // numazimuth
                        numAzim = i - (numScan * numazimuth)
                        Azimuth = buffer['azimuthal_axis'][numAzim]
                        self.dnames.append(f'#{numScan} Azm={Azimuth} {os.path.split(filename)[1][:60]}')
                    return True
            # test for MaxIV NeXus combined NXazint1d & NXazint2d
            test = True
            while test:  # block for combined NXazint1d and NXazint2d files,
                         # use break to bail out and try next block
                test = False
                entry = getNeXusBase(fp)
                subentry = getNeXusEntry(fp,entry,'NXsubentry')
                if len(subentry) == 0:
                    break # nothing to read
                for entry in subentry:
                    definition = fp[entry+'/definition'][()].decode()
                    if definition == 'NXazint1d' or definition == 'NXazint2d':
                        return True
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
        a buffer (if available) and subsequent calls can use the 
        buffer and will not need to access the file. 
        '''
        fpbuffer = kwarg.get('buffer',{})
        if not hasattr(self,'blknum'):
            if self.selections is None or len(self.selections) == 0:
                self.blknum = 0
            else:
                self.blknum = min(self.selections)
        # was file already read into buffer? If so, skip opening file to save time
        definition = fpbuffer.get('definition','')
        if definition == 'NXazint1d':
            return self.readNXazint1d(filename, fpbuffer)
        elif definition == 'NXazint2d':
            return self.readNXazint2d(filename, fpbuffer)

        try:        # first or non-buffered read
            fp = h5py.File(filename, 'r')
            entry = getNeXusBase(fp) # test for NeXus
            if entry:   # This is NeXus
                if 'definition' in fp[entry]: # MaxIV NXazint*
                    definition = fp[entry+'/definition'][()].decode()
                else: # is this a combined NXazint1d/NXazint2d file?
                    subentry = getNeXusEntry(fp,entry,'NXsubentry')
                    if len(subentry) == 0:
                        return False
                    elif len(subentry) == 1:
                        entry = subentry[0]
                    elif ParentFrame: # interactive, let the user decide
                        from .. import GSASIIctrlGUI as G2G
                        choices = ('NXazint1d 1D file','NXazint1d 2D file')
                        sel = G2G.ItemSelector(choices, ParentFrame=ParentFrame,
                                                   header='Select file section',
                                                   title='Select the section of the file to read')
                        if sel is None: return False
                        entry = subentry[sel]
                    else:   # scripted, assume if 2D is present, that is what is wanted
                        entry = subentry[1]
                    if 'definition' not in fp[entry]: return False
                    definition = fp[entry+'/definition'][()].decode()
                # got a file type, save it and if recognized, read it
                fpbuffer['definition'] = definition
                if definition == 'NXazint1d':
                    return self.readNXazint1d(filename, fpbuffer, entry)
                elif definition == 'NXazint2d':
                    return self.readNXazint2d(filename, fpbuffer, entry)
            return False # not a supported file type
        except IOError:  # unexpected since this was validated
            print (f'cannot open file {filename}')
            return False
        finally:
            fp.close()
        print (f'Unknown type of HDF5 powder file {filename}')
        return False

    def readNXazint1d(self, filename, fpbuffer={}, entry=None):
        '''Read HDF5 file in NeXus as produced by MAX IV as a NXazint1d.
        In this file, multiple scans are placed in a 2-D array (I and I_errors in 
        section NXdata), where one dimension is 2-theta and the other is a parametric 
        value such as temperature, time, etc. 

        see https://nxazint-hdf5-nexus-3229ecbd09ba8a773fbbd8beb72cace6216dfd5063e1.gitlab-pages.esrf.fr/classes/contributed_definitions/NXazint1d.html
        '''
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
        # test if we have what we need in the buffer and if not read it in
        if not self.readInNeXus(filename,fpbuffer,fileItems,'NXazint1d',entry): return False
        # now pull the selected dataset from the buffer
        self.numbanks = self.numparams 
        x = fpbuffer['radial_axis']
        y = fpbuffer['I'][self.blknum]
        try:
            esd = fpbuffer['I_errors'][self.blknum]
            w = np.where(esd==0,0,np.nan_to_num(1/esd**2))
        except:
            w = np.nan_to_num(1/y)    # best we can do, alas. W/o reported s.u.'s
        self.powderdata = [x,y,w,np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)]
        self.FillInParametics(fpbuffer,self.blknum)
        self.powderentry[0] = filename
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

    def readNXazint2d(self, filename, fpbuffer={}, entry=None):
        '''Read HDF5 file in NeXus as produced by MAX IV as a NXazint2d

        In this file, multiple scans are placed in a 3-D array (I and I_errors in 
        section NXdata), where one dimension is 2-theta and another is the azimuthal value
        and the third are a parametric value(s) such as temperature, time, etc.

        see https://nxazint-hdf5-nexus-3229ecbd09ba8a773fbbd8beb72cace6216dfd5063e1.gitlab-pages.esrf.fr/classes/contributed_definitions/NXazint2d.html
        '''
        self.comments = []
        fileItems = {
            # arrays
            'radial_axis':('NXdata','radial_axis'),
            'azimuthal_axis':('NXdata','azimuthal_axis'),
            'I':('NXdata','I'),
            'I_errors':('NXdata','I_errors'),
            # floats
            'wavelength':('NXmonochromator','wavelength'),
            'polarization_factor':('NXparameters','polarization_factor'),
            # strings
            'instrument/name':('NXinstrument','name'),
            'unit':('NXparameters','unit'),
            'azimuth_bins':('NXparameters','azimuth_bins'),
            'sample/name':('NXsample','name'),
            'source/name':('NXsource','name'),
        }
        # test if we have what we need in the buffer and if not read it in
        if not self.readInNeXus(filename,fpbuffer,fileItems,'NXazint2d',entry): return False
        # now pull the selected dataset from the buffer
        self.numazimuth = fpbuffer['azimuth_bins']
        self.numbanks = self.numparams * self.numazimuth
        # group by parametric variable
        numScan = self.blknum // self.numazimuth
        numAzim = self.blknum - (numScan * self.numazimuth)
        x = fpbuffer['radial_axis']
        y = fpbuffer['I'][numScan][numAzim]
        try:
            esd = fpbuffer['I_errors'][numScan][numAzim]
            w = np.where(esd==0,0,np.nan_to_num(1/esd**2))
        except:
            w = np.nan_to_num(1/y)    # best we can do, alas. W/o reported s.u.'s
        self.powderdata = [x,y,w,np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)]
        self.Sample['Azimuth'] = fpbuffer['azimuthal_axis'][numAzim]
        # add parametric var as a comment
        self.FillInParametics(fpbuffer,numScan)
        self.powderentry[0] = filename
        self.powderentry[2] = self.blknum  # bank number
        self.idstring = f'#{numScan} Azm={self.Sample["Azimuth"]} {os.path.split(filename)[1][:60]}'
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
    
    def readInNeXus(self,filename,fpbuffer,fileItems,fmt,entry=None):
        '''Read in items from NeXus labeled sections of the HDF5 file.

        For files where we are reading from a NXsubentry section
        rather than NXentry, variable `entry` is pointer to the 
        the selected NXsubentry section. If None, the NXentry
        is found. Otherwise `entry` points to the NXsubentry
        location, so only that portion of the tree is used.
        '''
        self.comments = []
        doread = False # has the file already been read into a buffer?
        for k in fileItems:
            if k not in fpbuffer:                
                doread = True
                break
        if doread:
            # Nope, need to fill the buffer
            try:
                fp = h5py.File(filename, 'r')
                if entry is None: entry = getNeXusBase(fp)
                # assemble list of used NeXus labels
                nexusDict = {i:None for i in set([i[0] for i in fileItems.values()])}
                # lookup keys for NeXus labels we will use
                recurseNeXusEntries(fp,entry,nexusDict)
                # save selected items from file into buffer
                # Convert all entries read into values or non-HDF5 objects so file
                # can be closed.
                savedKeys = [] # things that have already been read
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
                    elif 'int' in str(val.dtype):
                        fpbuffer[k] = int(val[()])
                    else:
                        fpbuffer[k] = val[()].decode()
                        self.comments.append(f'{k}={fpbuffer[k]}')
                if fpbuffer['unit'] != '2th':
                    print(f'{fmt} HDF5 file has units',fpbuffer['unit'])
                    self.errors = f'{fmt} only can be read with 2theta units'
                    return False
                self.numparams = len(fpbuffer['I'])
                # save arrays that are potentially tracking the parametric 
                # conditions into ParamTrackingVars. These arrays will have 
                # the same length as the number of datasets (self.numparams)
                if 'validate' not in fmt: # skip if we are validating the file rather than reading it
                    fpbuffer['ParamTrackingVars'] = {}
                    paramItems = []
                    for loc in nexusDict.values():
                        if loc is None: continue  # a NeXus label is not present
                        self.HDF5entries = []
                        paramItems = self.RecurseH5Element(fp[loc],length=self.numparams)
                        for i in paramItems:
                            for j in i:
                                key = loc+'/'+'/'.join(j)
                                if key in savedKeys: continue
                                savedKeys.append(key)
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
                print (f'Cannot open or read file {filename}')
                self.errors = f'{fmt} Can not open or read file {filename}'
                return False
            finally:
                fp.close()
            # initialize the block selection
            if self.selections is None or len(self.selections) == 0:
                self.blknum = 0
            else:
                self.blknum = min(self.selections)
        return True

    def FillInParametics(self,fpbuffer,count):
        '''put changing parametric variables into the comments
        '''
        for key,arr in fpbuffer['ParamTrackingVars'].items():
            val = arr[count]
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
    
    # HDF5 support routines.
    def RecurseH5Element(self,obj,prefix=[],length=None):
        '''Returns a list of entries of all keys in the HDF5 file
        (or group) in `obj`. Note that `obj` can be a file object, created by 
        `h5py.File` or can be a subsetgroup `fp['key/subkey']`.
        
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

    def ShowH5Element(self,obj,keylist):
        '''Format the contents of an HDF5 entry as a single line. Not used for 
        reading files, only used in :meth:`HDF5list`, which is here for software
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
    '''This returns a list of entries in a NeXus compilant HDF5 file matching 
    the name target, or an empty list, if this is not found. This only
    looks for the direct children of the key `base`.
    '''
    keyList = []
    for key in fp[base]:
        subkey = '/'.join([base,key])
        if "NX_class" in fp[subkey].attrs:
            #print(key, list(fp[subkey].attrs),fp[subkey].attrs["NX_class"])
            if ("NX_class" in fp[subkey].attrs and
                    fp[subkey].attrs["NX_class"] == target):
                keyList.append(subkey)
    return keyList

def recurseNeXusEntry(fp,node,target):
    '''Recurse through the HDF5 tree looking for NeXus class `target`.
    This stops after the first entry is found, and might be more useful
    if it returned a list when multiple definitions are present. 

    Not in use, as :func:`recurseNeXusEntries` is used to get all 
    targets in a single pass through the tree.
    '''
    if node is None: return  # needed?
    val = fp[node]
    if ("NX_class" in val.attrs and val.attrs["NX_class"] == target):
        return node
    if not isinstance(val, h5py.Group): return
    for key in val:
        subkey = '/'.join([node,key])
        res = recurseNeXusEntry(fp,subkey,target)
        if res: return res

def recurseNeXusEntries(fp,node,targetdict):
    '''recurse through the HDF5 tree looking for the NeXus classes
    in `targetdict`, storing the HDF5 key for each class in the dict.
    Note that if a NeXus classes is used more than once, only the 
    last encountered location will be saved. Use :func:`getNeXusEntry`
    when one needs to process all entries tagged with a class.

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
