# -*- coding: utf-8 -*-
'''
'''

from __future__ import division, print_function
import os
import sys
try:
    import h5py
except ImportError:
    h5py = None
import numpy as np
from .. import GSASIIobj as G2obj
from .. import GSASIIfiles as G2fil
#from .. import GSASIIpath

# things to do:
#   uncertainties
#   instr. parms
#instprmList = [('Bank',1.0), ('Lam',0.413263), ('Polariz.',0.99), 
#            ('SH/L',0.002), ('Type','PXC'), ('U',1.163), ('V',-0.126), 
#            ('W',0.063), ('X',0.0), ('Y',0.0), ('Z',0.0), ('Zero',0.0)]
#   comments
#   dataset naming
#   sample parameters
#sampleprmList = [('InstrName','APS 1-ID'), ('Temperature', 295.0)]
#  'Scale': [1.0, True], 'Type': 'Debye-Scherrer',
# 'Absorption': [0.0, False], 'DisplaceX': [0.0, False], 'DisplaceY': [0.0, False]# 'Pressure': 0.1, 'Time': 0.0, 'FreePrm1': 0.0,
# 'FreePrm2': 0.0, 'FreePrm3': 0.0, 'Gonio. radius': 200.0, 'Omega': 0.0,
# 'Chi': 0.0, 'Phi': 180.0, 'Azimuth': 0.0,
# 'Materials': [{'Name': 'vacuum', 'VolFrac': 1.0}, {'Name': 'vacuum', 'VolFrac': 0.0}],
# 'Thick': 1.0, 'Contrast': [0.0, 0.0], 'Trans': 1.0, 'SlitLen': 0.0}


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
            formatName = 'MAX IV HDF5',longFormatName = 'HDF5 integrated scans')
        self.scriptable = True
        #self.Iparm = {} #only filled for EDS data

    def ShowH5Element(self,obj,keylist):
        '''Format the contents of an HDF5 entry as a single line. Not used for 
        reading files, only used in :meth:`HDF5list`
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

    def RecurseH5Element(self,obj,prefix=[]):
        '''Returns a list of entries of all keys in the HDF5 file
        (or group) in `obj`. Note that `obj` can be a file object, created by 
        `h5py.File` or can be a subset `fp['key/subkey']`.

        The returned list is organized where: 
          * entry 0 is the top-level keys (/a, /b,...),
          * entry 1 has the first level keys (/a/c /a/d, /b/d, /b/e,...)
          * ...
        Not used for reading files, used only in :meth:`HDF5list`
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
            self.HDF5entries[depth].append(nextprefix)
            # check for link objects
            l = obj.get(i, getlink=True)
            if isinstance(l, h5py.ExternalLink): continue
            try:
                typ = str(type(obj[i]))
            except:
                print(f'**Error** with key {prefix}/{i}')
                continue
            if ".Group'" in typ:
                #t = f'{prefix}/{i}'
                #print(f'\n{nextprefix} contents {(60-len(t))*'='}')
                self.RecurseH5Element(obj[i],nextprefix)
        return self.HDF5entries
        
                
    def HDF5list(self, filename):
        '''Shows the contents of an HDF5 file as a short listing. 
        This is not used for HDF5 reading, but is of help with a new
        type of HDF5 file to see what is present.

        :param filename: 
        '''
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
                lbl = self.ShowH5Element(fp,k)
                if '\n' in lbl:
                    lbl = '; '.join(lbl.split('\n'))
                if len(lbl) > 50:
                    lbl = lbl[:50] + '...'
                # if '\n' in lbl:
                #     lbl = lbl.split()[0] + '...'
                if lbl != '(group)': strings.append(f"{'/'.join(k):{m}s} {lbl}")
        with open(filename+'_contents.txt', 'w') as fp:
            for i in strings: fp.write(f'{i}\n')
                    
    def ContentsValidator(self, filename):
        '''Test if valid by seeing if the HDF5 library recognizes the file. 
        Then get file type (currently MAX IV NeXus/NXazint[12]d only)
        '''
        #from .. import GSASIIpath
        try:
            fp = h5py.File(filename, 'r')
            if 'entry' in fp: # NeXus
                #self.HDF5entries = []
                #self.HDF5list(filename)
                if 'definition' in fp['/entry']:
                    # MAX IV NXazint1d file
                    if fp['/entry/definition'][()].decode() == 'NXazint1d':
                        return True
                    # MAX IV NXazint1d file
                    #if fp['/entry/definition'][()].decode() == 'NXazint2d':
                    #    return True
        except IOError:
            return False
        finally:
            fp.close()
        return False

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
        try:
            fp = h5py.File(filename, 'r')
            if 'entry' in fp: # NeXus
                if 'definition' in fp['/entry']:
                    # MAX IV NXazint1d file
                    if fp['/entry/definition'][()].decode() == 'NXazint1d':
                        return self.readNXazint1d(filename, fpbuffer)

                    # MAX IV NXazint1d file
                    #if fp['/entry/definition'][()].decode() == 'NXazint2d':
                    #    return self.readNXazint2d(filename, fpbuffer)
                    #    return True
                    # https://nxazint-hdf5-nexus-3229ecbd09ba8a773fbbd8beb72cace6216dfd5063e1.gitlab-pages.esrf.fr/classes/contributed_definitions/NXazint2d.html
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
        doread = False # has the file already been read into a buffer?
        arrays = ('entry/data/radial_axis','entry/data/I')
        floats = ('entry/instrument/monochromator/wavelength',
                  'entry/reduction/input/polarization_factor')
        strings = ('entry/instrument/source/name','entry/reduction/input/unit')
        for i in arrays+floats+strings:
            if i not in fpbuffer:
                doread = True
                break
        if doread:   # read into buffer
            try:
                fp = h5py.File(filename, 'r')
                for i in arrays:
                    fpbuffer[i] = np.array(fp.get(i))
                for i in floats:
                    fpbuffer[i] = float(fp[i][()])
                for i in strings:
                    fpbuffer[i] = fp[i][()].decode()
                if fpbuffer['entry/reduction/input/unit'] != '2th':
                    print('NXazint1d HDF5 file has units',fpbuffer['entry/reduction/input/unit'])
                    self.errors = 'NXazint1d only can be read with 2th units'
                    return False
                if self.selections is None or len(self.selections) == 0:
                    self.blknum = 0
                else:
                    self.blknum = min(self.selections)
            except IOError:
                print ('cannot open file '+ filename)
                return False
            finally:
                fp.close()
            self.numbanks=len(fpbuffer['entry/data/I'])
            # # get overriding sample & instrument parameters
            # fpbuffer['sampleprm'] = {}
            # samplefile = os.path.splitext(filename)[0] + '.samprm'
            # if os.path.exists(samplefile):
            #     fp = open(samplefile,'r')
            #     S = fp.readline()
            #     while S:
            #         if not S.strip().startswith('#'):
            #             [item,val] = S[:-1].split(':')
            #             fpbuffer['sampleprm'][item.strip("'")] = eval(val)
            #         S = fp.readline()
            #     fp.close()
            # fpbuffer['instprm'] = {}
            # instfile = os.path.splitext(filename)[0] + '.instprm'
            # if os.path.exists(instfile):
            #     self.instmsg = 'HDF and .instprm files'
            #     fp = open(instfile,'r')
            #     S = fp.readline()
            #     while S:
            #         if not S.strip().startswith('#'):
            #             [item,val] = S[:-1].split(':')
            #             fpbuffer['instprm'][item.strip("'")] = eval(val)
            #         S = fp.readline()
            #     fp.close()
        # now transfer information into current histogram 
        #self.pwdparms['Instrument Parameters'] = [
        #    {'Type': ['PXC', 'PXC', False]},
        #    {}]
        # inst = {}
        # inst.update(instprmList)
        # inst.update(fpbuffer['instprm'])
        # for key,val in inst.items():
        #     self.pwdparms['Instrument Parameters'][0][key] = [val,val,False]
        # samp = {}
        # samp.update(sampleprmList)
        # samp.update(fpbuffer['sampleprm'])
        # for key,val in samp.items():
        #     self.Sample[key] = val
        x = fpbuffer['entry/data/radial_axis']
        y = fpbuffer['entry/data/I'][self.blknum]
        w = np.nan_to_num(1/y)    # this is not correct
        #self.pwdparms['Instrument Parameters'][0]['Azimuth'] = [90-eta,90-eta,False]
        #self.pwdparms['Instrument Parameters'][0]['Bank'] = [self.blknum,self.blknum,False]
#        self.Sample['Gonio. radius'] = float(S.split('=')[1])
#        self.Sample['Omega'] = float(S.split('=')[1])
#        self.Sample['Chi'] = float(S.split('=')[1])
        #self.Sample['Phi'] = Omega = fpbuffer['Omegas'][self.blknum]
        self.powderdata = [x,y,w,np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)]
        #self.comments = comments[selblk]
        self.powderentry[0] = filename
        #self.powderentry[1] = Pos # position offset (never used, I hope)
        self.powderentry[2] = self.blknum  # bank number
        self.idstring = f'#{self.blknum} {os.path.split(filename)[1][:60]}'
        self.instdict['wave'] = fpbuffer['entry/instrument/monochromator/wavelength']
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
