# -*- coding: utf-8 -*-
'''Import a collection of "lineouts" from MIDAS from a zarr zip file
'''

from __future__ import division, print_function
import os,sys
try:
    import zarr
except ImportError:
    zarr = None
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath

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
    midassections = ('InstrumentParameters','Omegas', 'REtaMap', 'OmegaSumFrame')
    def __init__(self):
        if zarr is None:
            self.UseReader = False
            print('MIDAS_Zarr Reader skipped because zarr module is not installed')
            os.path.split(sys.executable)[0]
            conda = os.path.join(os.path.split(sys.executable)[0],'conda')
            if os.path.exists(conda):
                print(f'To fix this use command:\n\t{conda} install zarr')
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.zarr.zip',),strictExtension=True,
            formatName = 'MIDAS zarr',longFormatName = 'MIDAS zarr intergrated scans')
        self.scriptable = True
        #self.Iparm = {} #only filled for EDS data

    def ContentsValidator(self, filename):
        '''Test if valid by seeing if the zarr module recognizes the file. Then
        get file type (currently Midas only)
        '''
        try:
            fp = zarr.open(filename, 'r')
            if all([(i in fp) for i in self.midassections]): # are expected MIDAS sections present?
                self.mode = 'midas'
                return True # must be present for Midas output
        except:
            return False
        finally:
            del fp
        return False

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
        # TODO: need to think about normalization
        # sample parameters (such as temperature)
        # xfer instrument from fpbuffer
        # review xfer of sample parameters
        # 
        fpbuffer = kwarg.get('buffer',{})
        if not hasattr(self,'blknum'):
            self.blknum = 0    # image counter for multi-image files
        # check if this is a valid MIDAS file
        if self.mode is None:
            try:
                fp = zarr.open(filename, 'r')
                 # are expected MIDAS sections present?
                if all([(i in fp) for i in self.midassections]):
                    self.mode = 'midas'
                else:
                    print (f'{filename} is not a MIDAS file')
                    return False
            except:
                print (f'cannot read as zarr file: {filename}')
                return False
            finally:
                del fp
                    
        if self.mode == 'midas':
            return self.readMidas(filename, fpbuffer)

        return False

    def readMidas(self, filename, fpbuffer={}):
        '''Read zarr file produced by Midas
        '''
        self.instmsg = 'MIDAS zarr file'
        # has the zarr file already been cached?
        doread = False # yes
        for i in ('intenArr','REtaMap','omegas', 'REtaMap', 'unmasked', '2Read'):
            if i not in fpbuffer:
                doread = True # no
                break

        #======================================================================
        # cache the contents of the zarr file on the first call to this
        # (or every call if no buffer is supplied -- very slow)
        #======================================================================
        if doread:   # read into buffer
            print()
            try:
                fp = zarr.open(filename, 'r')
                fpbuffer['REtaMap'] = np.array(fp['REtaMap']) # 4 x Nbins x Nazimuth
                # [0]: radius; [1] 2theta; [2] eta; [3] bin area 
                # tabulate integrated intensities image & eta
                fpbuffer['intenArr'] = []
                fpbuffer['omegas'] = []
                for i,k in enumerate(fp['OmegaSumFrame']):
                    fpbuffer['intenArr'].append(fp['OmegaSumFrame'][k])
                    fpbuffer['omegas'].append(0.5*
                        (fp['OmegaSumFrame']['LastFrameNumber_1'].attrs['FirstOme']+
                         fp['OmegaSumFrame']['LastFrameNumber_1'].attrs['LastOme']))                    
                Nimg = len(fp['OmegaSumFrame'])   # number of images
                Nbins,Nazim = fpbuffer['REtaMap'][1].shape
                # Nbins: number of points in each lineout (not all of which may be
                #     used due to geometrical masking)
                # Nazim: number of azimuthal "cake slices"

                # get a list of points in use at each azimuth
                unmasked = fpbuffer['unmasked'] = [(fpbuffer['REtaMap'][3][:,i] != 0) for i in range(Nazim)]
                # find the azimuths with more than 20 points
                mAzm = [i for i in range(Nazim) if sum(unmasked[i]) > 20]
                
                # generate a list of lineouts to be read from selections
                if self.selections:
                    sel = self.selections
                else:
                    sel = range(Nimg*Nazim) # nothing selected, use all
                fpbuffer['2Read'] = [(i%Nazim,i//Nazim) for i in sel if i%Nazim in mAzm]
                # fpbuffer['2Read'] is the list of lineouts to be read, where each entry
                # contains two values, the azumuth and the image number (iAzm,iImg)
                # defined points for each lineout are then 
                #   intensities : fpbuffer['intenArr'][iImg][:,iAzm][unmasked[iAzm]]
                #   2thetas: fpbuffer['REtaMap'][1][:,iAzm][unmasked[iAzm]]
                self.MIDASinstprm = {i:j[0] # reform as a native dict
                    for i,j in fp['InstrumentParameters'].items()} 
            except IOError:
                print ('cannot open file '+ filename)
                return False
            finally:
                del fp
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
            if os.path.exists(instfile):
                self.instmsg = 'zarr and .instprm files'
                fp = open(instfile,'r')
                S = fp.readline()
                while S:
                    if not S.strip().startswith('#'):
                        [item,val] = S[:-1].split(':')
                        fpbuffer['instprm'][item.strip("'")] = eval(val)
                    S = fp.readline()
                del fp
        #======================================================================
        # start reading 
        #======================================================================
        # look for the next non-empty scan (lineout)
        iAzm,iImg = fpbuffer['2Read'][self.blknum]
        unmasked = fpbuffer['unmasked']
        y = fpbuffer['intenArr'][iImg][:,iAzm][unmasked[iAzm]]
        x = fpbuffer['REtaMap'][1][:,iAzm][unmasked[iAzm]]
        normalization = fpbuffer['REtaMap'][3][:,iAzm][unmasked[iAzm]]
        # note that y values have been scaled by division by normalization values
        # esd = np.sqrt(y*normalization)
        w = np.nan_to_num(1./(y * normalization))
        omega = fpbuffer['omegas'][iImg]
        eta = sum(fpbuffer['REtaMap'][2][:,iAzm][unmasked[iAzm]])/sum(
            unmasked[iAzm]) # single-valued?
        radius = 1000
#        radius = sum(fpbuffer['REtaMap'][0][:,iAzm][unmasked[iAzm]])/sum(
#            unmasked[iAzm]) * self.pixelsize
        
        # now transfer instprm & sample prms into current histogram 
        self.pwdparms['Instrument Parameters'] = [{}, {}]
        inst = {}
        inst.update(instprmList)
        inst.update(self.MIDASinstprm)
        for key,val in inst.items():
            self.pwdparms['Instrument Parameters'][0][key] = [val,val,False]
        samp = {}
        samp.update(sampleprmList)
        samp.update(self.MIDASsampleprm)
        for key,val in samp.items():
            self.Sample[key] = val
        self.numbanks=len(fpbuffer['2Read'])  # number of lineouts to be read
        self.pwdparms['Instrument Parameters'][0]['Azimuth'] = [90-eta,90-eta,False]
        self.pwdparms['Instrument Parameters'][0]['Bank'] = [iAzm,iAzm,False]
        self.Sample['Gonio. radius'] = float(radius)
#        self.Sample['Omega'] = float(S.split('=')[1])
#        self.Sample['Chi'] = float(S.split('=')[1])
        self.Sample['Phi'] = omega
#        self.Sample['FreePrm1']
#        self.Controls['FreePrm1'] = 'Midas Label'  # relabel the parameter
#        self.Sample['Temperature']
#        self.Sample['Pressure']
#        self.Sample['Time']
        self.powderdata = [x,y,w,np.zeros_like(x),np.zeros_like(x),np.zeros_like(x)]
        #self.comments = comments[selblk]
        self.powderentry[0] = filename
        #self.powderentry[1] = Pos # position offset (never used, I hope)
        self.powderentry[2] = self.blknum  # consecutive bank number
        self.idstring = f'{os.path.split(filename)[1][:10]} img={iImg} omega={omega} eta={eta}'
        if GSASIIpath.GetConfigValue('debug'):
            print(f'Read entry #{iAzm} img# {iImg} from file {filename}\n{self.idstring}')
        else:
            print('.',end='')  # do something since this can take a while
            sys.stdout.flush()
        
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
        return True
