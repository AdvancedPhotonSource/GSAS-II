# -*- coding: utf-8 -*-
'''Import a collection of "lineouts" from MIDAS from a zarr zip file
'''
# TODO: radius to be added to Zarr file
#
from __future__ import division, print_function
import os
try:
    import zarr
except ImportError:
    zarr = None
import numpy as np
from .. import GSASIIobj as G2obj
from .. import GSASIIfiles as G2fil
from .. import GSASIIpath

instprmList = [('Bank',1.0), ('Lam',0.413263), ('Polariz.',0.99),
            ('SH/L',0.002), ('Type','PXC'), ('U',1.163), ('V',-0.126),
            ('W',0.063), ('X',0.0), ('Y',0.0), ('Z',0.0), ('Zero',0.0)]
#   comments
#   sample parameters
sampleprmList = [('InstrName','APS 1-ID'), ('Temperature', 295.0)]
#  'Scale': [1.0, True], 'Type': 'Debye-Scherrer',
# 'Absorption': [0.0, False], 'DisplaceX': [0.0, False], 'DisplaceY': [0.0, False]# 'Pressure': 0.1, 'Time': 0.0, 'FreePrm1': 0.0,
# 'FreePrm2': 0.0, 'FreePrm3': 0.0, 'Gonio. radius': 200.0, 'Omega': 0.0,
# 'Chi': 0.0, 'Phi': 180.0, 'Azimuth': 0.0,
# 'Materials': [{'Name': 'vacuum', 'VolFrac': 1.0}, {'Name': 'vacuum', 'VolFrac': 0.0}],
# 'Thick': 1.0, 'Contrast': [0.0, 0.0], 'Trans': 1.0, 'SlitLen': 0.0}

class MIDAS_Zarr_Reader(G2obj.ImportPowderData):
    '''Routine to read multiple powder patterns from a Zarr file
    created by MIDAS. Perhaps in the future, other software might also
    use this file type as well.

    For Midas, the main file is <file>.zip, but optionally sample and
    instrument parameters can be placed in <file>.samprm and <file>.instprm.
    Any parameters placed in those files will override values set in the .zip
    file.
    '''
    mode = None
    midassections = ('InstrumentParameters', 'REtaMap', 'OmegaSumFrame')
    def __init__(self):
        if zarr is None:
            self.UseReader = False
            msg = 'MIDAS_Zarr Reader skipped because zarr module is not installed.'
            G2fil.ImportErrorMsg(msg,{'MIDAS Zarr importer':['zarr']})
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.zarr.zip',),strictExtension=True,
            formatName = 'MIDAS zarr',longFormatName = 'MIDAS zarr intergrated scans')
        self.scriptable = True
        #self.Iparm = {} #only filled for EDS data

    def zarrOpen(self,filename):
        import asyncio
        async def zarrAsyncOpen(filename):
            fp = store = None
            try:
                fp = zarr.open(filename, mode='r')
                return fp,store
            except:
                # workaround for zarr 3.0.x where "auto discover" is
                # not yet implemented
                # (https://github.com/zarr-developers/zarr-python/issues/2922)
                try:
                    store = await zarr.storage.ZipStore.open(filename, mode="r")
                    fp = zarr.open_group(store, mode='r')
                    return fp,store
                except FileNotFoundError as msg:
                    print (f'cannot read as zarr file: {filename}')
                except Exception as msg:
                    self.errors = f'Exception from zarr module (version={zarr.__version__}):'
                    self.errors += '\n\t' + str(msg)
            return None,None
        fp,store = asyncio.run(zarrAsyncOpen(filename))
        return fp,store
                        
    def ContentsValidator(self, filename):
        '''Test if valid by seeing if the zarr module recognizes the file. Then
        get file type (currently Midas only)
        '''
        fp, store = self.zarrOpen(filename)
        if not fp:
            self.errors = f'Exception from zarr module (version={zarr.__version__}):'
            self.errors += '\n\t' + str(msg)
            self.errors += '\nPlease report this'
        elif all([(i in fp) for i in self.midassections]): # are expected MIDAS sections present?
            self.mode = 'midas'
            del fp, store
            return True # Passes test for Midas output
        #else:  # test here for any other formats that might use zarr
        #    pass
        del fp, store
        return

    def Reader(self, filename, ParentFrame=None, **kwarg):
        '''Scan file for sections needed by defined file types (currently
        only Midas) and then use appropriate routine to read the file.
        Most of the time the results are placed in the buffer arg (if supplied)
        so the file is not read on most calls.

        For MIDAS, masking can eliminate some or all points in an azimuthal
        region. This will only return "lineouts" (aka diffractograms) that
        have 20 or more points in them.

        Note that if Reader.selections is used to select individual
        "lineouts", the selections are numbered against all possible
        "lineouts" not the ones that have 20 or more points.
        '''
        fpbuffer = kwarg.get('buffer',{})
        if not hasattr(self,'blknum'):
            self.blknum = 0    # image counter for multi-image files
        # check if this is a valid MIDAS file
        if self.mode is None:
            fp, store = self.zarrOpen(filename)
            # are expected MIDAS sections present?
            if all([(i in fp) for i in self.midassections]):
                self.mode = 'midas'
            else:
                print (f'{filename} is not a MIDAS file')
            del fp, store
        #else:  # test here for any other formats that might use zarr
        #    pass
        
        if self.mode == 'midas':
            res = self.readMidas(filename, fpbuffer)
        else:
            res = False
        return res

    def readMidas(self, filename, fpbuffer={}):
        '''Read zarr file produced by Midas
        '''
        self.instmsg = 'MIDAS zarr file'
        # has the zarr file already been cached?
        doread = False # yes
        for i in ('intenArr','REtaMap','attributes', 'REtaMap', 'unmasked', '2Read'):
            if i not in fpbuffer:
                doread = True # no
                break

        #======================================================================
        # cache the contents of the zarr file on the first call to this
        # (or every call if no buffer is supplied -- very slow)
        #======================================================================
        if doread:   # read into buffer
            print('Caching MIDAS zarr contents...')
            try:
                fp, store = self.zarrOpen(filename)
                fpbuffer['REtaMap'] = np.array(fp['REtaMap']) # 4 x Nbins x Nazimuth
                # [0]: radius; [1] 2theta; [2] eta; [3] bin area
                # tabulate integrated intensities image & eta
                fpbuffer['intenArr'] = []
                fpbuffer['attributes'] = []
                for i,k in enumerate(fp['OmegaSumFrame']):
                    fpbuffer['intenArr'].append(fp['OmegaSumFrame'][k])
                    fpbuffer['attributes'].append(
                        dict(fp['OmegaSumFrame'][k].attrs.items()))
                Nimg = len(fp['OmegaSumFrame'])   # number of images
                Nbins,Nazim = fpbuffer['REtaMap'][1].shape
                # Nbins: number of points in each lineout (not all of which may be
                #     used due to geometrical masking)
                # Nazim: number of azimuthal "cake slices"

                # get a list of points in use at each azimuth
                fpbuffer['unmasked'] = [(fpbuffer['REtaMap'][3][:,i] != 0) for i in range(Nazim)] # will be True if area is >0
                # find the azimuths with more than 20 points
                mAzm = [i for i in range(Nazim) if sum(fpbuffer['unmasked'][i]) > 20]

                # generate a list of lineouts to be read from selections
                if self.selections:
                    sel = self.selections
                else:
                    sel = range(Nimg*Nazim) # nothing selected, use all
                # fpbuffer['2Read'] is the list of lineouts to be read, where each entry
                # contains two values, the azumuth and the image number (iAzm,iImg)
                # defined points for each lineout are then
                #   intensities : fpbuffer['intenArr'][iImg][:,iAzm][unmasked[iAzm]]
                #   2thetas: fpbuffer['REtaMap'][1][:,iAzm][unmasked[iAzm]]
                fpbuffer['2Read'] = [(i%Nazim,i//Nazim) for i in sel if i%Nazim in mAzm]
                # xfrom Zarr dict into a native dict
                self.MIDASinstprm = {i:fp['InstrumentParameters'][i][0] for i in fp['InstrumentParameters']}
                # change a few keys
                for key,newkey in [('Polariz','Polariz.'),('SH_L','SH/L')]:
                    if key in self.MIDASinstprm:
                        self.MIDASinstprm[newkey] = self.MIDASinstprm[key]
                        del self.MIDASinstprm[key]
            except IOError:
                print ('cannot open file '+ filename)
                return False
            finally:
                del fp, store
            # get overriding sample & instrument parameters
            self.MIDASsampleprm = {}
            samplefile = os.path.splitext(filename)[0] + '.samprm'
            if os.path.exists(samplefile):
                fp = open(samplefile,'r')
                S = fp.readline()
                while S:
                    if not S.strip().startswith('#'):
                        [item,val] = S[:-1].split(':')
                        self.MIDASsampleprm[item.strip("'")] = eval(val)
                    S = fp.readline()
                fp.close()
            # overload from file, if found
            instfile = os.path.splitext(filename)[0] + '.instprm'
            self.instvals = [{},{}]
            if os.path.exists(instfile):
                self.instmsg = 'zarr and .instprm files'
                fp = open(instfile,'r')
                instLines = fp.readlines()
                fp.close()
                nbank,self.instvals = G2fil.ReadInstprm(instLines, None, self.Sample)
        #======================================================================
        # start reading
        #======================================================================
        # look for the next non-empty scan (lineout)
        iAzm,iImg = fpbuffer['2Read'][self.blknum]
        nFrame = fpbuffer['attributes'][iImg].get('Number Of Frames Summed',1.)
        unmasked = fpbuffer['unmasked'][iAzm]
        y = fpbuffer['intenArr'][iImg][:,iAzm][unmasked]/nFrame
        x = fpbuffer['REtaMap'][1][:,iAzm][unmasked]
        normalization = fpbuffer['REtaMap'][3][:,iAzm][unmasked]
        # compute the uncertainty on the normalized intensities
        #   Y(normalized) = Y-norm = Ysum / area
        #   sigma(Y-norm) = sqrt(Ysum) / area = sqrt[Y-norm * area] / area
        #                 = sqrt( Y-norm / area)
        # For GSAS-II want to normalize by the number of frames,
        #   Y(GSAS) = Y-norm / nFrame
        #   sigma(Y-GSAS) = sigma(Y-norm) / nFrame
        #                 = sqrt( Y-norm / area) / nFrame
        #   weight(Y-GSAS) is 1/sigma[Y-GSAS]**2 = nFrame**2 * area / Y-norm
        w = np.where(y > 0, np.zeros_like(y), nFrame**2 * normalization/ y )
        omega = 999.  # indicates an error
        try:
            omega = 0.5 * (
                fpbuffer['attributes'][iImg]['FirstOme'] +
                fpbuffer['attributes'][iImg]['LastOme'])
        except:
            pass
        eta = sum(fpbuffer['REtaMap'][2][:,iAzm][unmasked])/sum(unmasked) # compute an averaged Phi value
        radius = 1000   # sample to detector distance (mm)
        if 'Distance' in self.MIDASinstprm:
            radius = self.MIDASinstprm['Distance']/1000   # convert from microns

        # now transfer instprm & sample prms into current histogram
        self.pwdparms['Instrument Parameters'] = [{}, {}]
        inst = {}
        inst.update(instprmList)
        for i in inst:
            if i in self.MIDASinstprm:
                inst[i] = self.MIDASinstprm[i]
        for key,val in inst.items():
            self.pwdparms['Instrument Parameters'][0][key] = [val,val,False]
        self.pwdparms['Instrument Parameters'][0].update(self.instvals[0])
        self.pwdparms['Instrument Parameters'][0]['Azimuth'] = [90-eta,90-eta,False]
        self.pwdparms['Instrument Parameters'][0]['Bank'] = [iAzm,iAzm,False]
        samp = {}
        samp.update(sampleprmList)
        samp.update(self.MIDASsampleprm)
        try:
            sampleprmList['Temperature'] = fpbuffer['attributes'][iImg]['Temperature']
        except:
            pass
        try:
            sampleprmList['Pressure'] = fpbuffer['attributes'][iImg]['Pressure']
        except:
            pass
        for key,val in samp.items():
            self.Sample[key] = val
        self.numbanks=len(fpbuffer['2Read'])  # number of remaining lineouts to be read
        # save the various GSAS-II instrument and sample parameters
        self.Sample['Gonio. radius'] = float(radius)
#        self.Sample['Omega'] = float(S.split('=')[1])
#        self.Sample['Chi'] = float(S.split('=')[1])
        self.Sample['Phi'] = omega
#        self.Sample['FreePrm1']
#        self.Controls['FreePrm1'] = 'Midas Label'  # relabel the parameter
#        self.Sample['Time']
        self.powderdata = [x,y,w,np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)]
        #self.comments = comments[selblk]
        self.powderentry[0] = filename
        #self.powderentry[1] = Pos # position offset (never used, I hope)
        self.powderentry[2] = self.blknum  # consecutive bank number
        self.idstring = f'{os.path.split(filename)[1]} img={iImg} omega={omega} eta={eta}'
        if GSASIIpath.GetConfigValue('debug'):
            print(f'Read entry #{iAzm} img# {iImg} from file {filename}\n{self.idstring}')
#        else:
#            print('.',end='')  # do something since this can take a while
#            sys.stdout.flush()
        # are there more lineouts after this one in current image to read?
        self.blknum += 1
        if self.blknum >= len(fpbuffer['2Read']): # is this the last scan?
            self.repeat = False
        else:
            self.repeat = True
            return True
        # read last pattern, reset buffer & counter
        print()
        self.blknum = 0
        for key in list(fpbuffer.keys()): del fpbuffer[key]
        print('...read done')
        return True
