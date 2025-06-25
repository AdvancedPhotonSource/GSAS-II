# -*- coding: utf-8 -*-
'''Class to read a phase from a CIF
'''
from __future__ import division, print_function
import numpy as np
import os.path
from .. import GSASIIpath
from .. import GSASIIobj as G2obj
from .. import GSASIIfiles as G2fil
try:
    import CifFile as cif # PyCifRW from James Hester as a package
except ImportError:
    try:
        from .. import CifFile as cif # PyCifRW, as distributed w/G2 (old)
    except ImportError:
        cif = None
asind = lambda x: 180.*np.arcsin(x)/np.pi

class CIFpwdReader(G2obj.ImportPowderData):
    'Routines to import powder data from a CIF file'
    def __init__(self):
        if cif is None:
            self.UseReader = False
            msg = 'CIF PWDR Reader skipped because PyCifRW (CifFile) module is not installed.'
            G2fil.ImportErrorMsg(msg,{'CIF powder importer':['pycifrw']})
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.CIF','.cif'),
            strictExtension=False,
            formatName = 'CIF',
            longFormatName = 'Powder data from CIF'
            )
    # Validate the contents
    def ContentsValidator(self, filename):
        'Use standard CIF validator'
        fp = open(filename,'r')
        return self.CIFValidator(fp)
        fp.close()

    def Reader(self,filename, ParentFrame=None, **kwarg):
        '''Read powder data from a CIF.
        If multiple datasets are requested, use self.repeat and buffer caching.
        '''

        # Define lists of data names used for holding powder diffraction data
        # entries of a type that are not implemented are commented out in case
        # we will want them later.
        xDataItems = (  # "x-axis" data names"
            ('_pd_meas_2theta_range_min', '_pd_meas_2theta_range_max', '_pd_meas_2theta_range_inc'),
            ('_pd_proc_2theta_range_min', '_pd_proc_2theta_range_max', '_pd_proc_2theta_range_inc'),
            '_pd_meas_2theta_scan',
            '_pd_meas_time_of_flight',
            '_pd_proc_2theta_corrected',
            #'_pd_proc_d_spacing',
            #'_pd_proc_energy_incident',
            #'_pd_proc_energy_detection',
            '_pd_proc_recip_len_q',
            #'_pd_proc_wavelength',
        )
        intDataItems = ( # intensity axis data names
            '_pd_meas_counts_total',
            #'_pd_meas_counts_background',
            #'_pd_meas_counts_container',
            '_pd_meas_intensity_total',
            #'_pd_meas_intensity_background',
            #'_pd_meas_intensity_container',
            '_pd_proc_intensity_net',
            '_pd_proc_intensity_total',
            #'_pd_proc_intensity_bkg_calc',
            #'_pd_proc_intensity_bkg_fix',
            '_pd_calc_intensity_net',   # allow computed patterns as input data?
            '_pd_calc_intensity_total',
            )

        ESDDataItems = ( # items that can be used to compute uncertainties for obs data
            '_pd_proc_ls_weight',
            '_pd_meas_counts_total'
            )

        ModDataItems = ( # items that modify the use of the data
            '_pd_meas_step_count_time',
            '_pd_meas_counts_monitor',
            '_pd_meas_intensity_monitor',
            '_pd_proc_intensity_norm',
            '_pd_proc_intensity_incident',
            )
        rdbuffer = kwarg.get('buffer')
        cf = None
        choicelist = None
        selections = None
        # reload previously saved values
        if self.repeat and rdbuffer is not None:
            cf = rdbuffer.get('lastcif')
            choicelist = rdbuffer.get('choicelist')
            print ('debug: Reuse previously parsed CIF')
            selections = rdbuffer.get('selections')
        if cf is None:
            if GSASIIpath.GetConfigValue('debug'): print("Starting parse of {} as CIF".format(filename))
            cf = G2obj.ReadCIF(filename)
            if GSASIIpath.GetConfigValue('debug'): print ("CIF file parsed")
        # scan all blocks for sets of data
        if choicelist is None:
            choicelist = []
            for blk in cf.keys():
                blkkeys = [k.lower() for k in cf[blk].keys()] # save a list of the data items, since we will use it often
                # scan through block for x items
                xldict = {}
                for x in xDataItems:
                    if type(x) is tuple: # check for the presence of all three items that define a range of data
                        if not all([i in blkkeys for i in x]): continue
                        try:
                            items = [float(cf[blk][xi]) for xi in x]
                            l = 1 + int(0.5 + (items[1]-items[0])/items[2])
                        except:
                            continue
                    else:
                        if x not in blkkeys: continue
                        l = len(cf[blk][x])
                    if xldict.get(l) is None:
                        xldict[l] = [x]
                    else:
                        xldict[l].append(x)
                # now look for matching intensity items
                yldict = {}
                suldict = {}
                for y in intDataItems:
                    if y in blkkeys:
                        l = len(cf[blk][y])
                        if yldict.get(l) is None:
                            yldict[l] = [y]
                        else:
                            yldict[l].append(y)
                        # now check if the first item has an uncertainty
                        if cif.get_number_with_esd(cf[blk][y][0])[1] is None: continue
                        if suldict.get(l) is None:
                            suldict[l] = [y]
                        else:
                            suldict[l].append(y)
                for y in ESDDataItems:
                    if y in blkkeys:
                        l = len(cf[blk][y])
                        if suldict.get(l) is None:
                            suldict[l] = [y]
                        else:
                            suldict[l].append(y)
                modldict = {}
                for y in ModDataItems:
                    if y in blkkeys:
                        l = len(cf[blk][y])
                        if modldict.get(l) is None:
                            modldict[l] = [y]
                        else:
                            modldict[l].append(y)
                for l in xldict:
                    if yldict.get(l) is None: continue
                    choicelist.append([blk,l,xldict[l],yldict[l],suldict.get(l,[]),modldict.get(l,[])])
                    #print blk,l,xldict[l],yldict[l],suldict.get(l,[]),modldict.get(l,[])
            print ("CIF file scanned for blocks with data")
        if not choicelist:
            selblk = None # no block to choose
            self.errors = "No powder diffraction blocks found"
            return False
        elif len(choicelist) == 1: # only one choice
            selblk = 0
        elif self.repeat and selections is not None:
            # we were called to repeat the read
            print ('debug: repeat #',self.repeatcount,'selection',selections[self.repeatcount])
            selblk = selections[self.repeatcount]
            self.repeatcount += 1
            if self.repeatcount >= len(selections): self.repeat = False
        else:                       # choose from options
            # compile a list of choices for the user
            choices = []
            for blk,l,x,y,su,mod in choicelist:
                sx = x[0]
                if len(x) > 1: sx += '...'
                sy = y[0]
                if len(y) > 1: sy += '...'
                choices.append(
                    'Block '+str(blk)+', '+str(l)+' points. X='+sx+' & Y='+sy
                    )
            from .. import GSASIIctrlGUI as G2G
            selections = G2G.MultipleBlockSelector(
                choices,
                ParentFrame=ParentFrame,
                title='Select dataset(s) to read from the list below',
                size=(600,100),
                header='Dataset Selector')
            if len(selections) == 0:
                self.errors = "Abort: block not selected"
                return False
            selblk = selections[0] # select first in list
            if len(selections) > 1: # prepare to loop through again
                self.repeat = True
                self.repeatcount = 1
                if rdbuffer is not None:
                    rdbuffer['selections'] = selections
                    rdbuffer['lastcif'] = cf # save the parsed cif for the next loop
                    rdbuffer['choicelist'] = choicelist # save the parsed choices for the future

        # got a selection, now read it
        # do we need to ask which fields to read?
        blk,l,xch,ych,such,modch = choicelist[selblk]
        xi,yi,sui,modi = 0,0,0,0
        if len(xch) > 1 or len(ych) > 1 or len(such) > 1 or len(modch) > 0:
            choices = []
            chlbls = []
            chlbls.append('Select the scanned data item')
            xchlst = []
            for x in xch:
                if type(x) is tuple:
                    xchlst.append(x[0])
                else:
                    xchlst.append(x)
            choices.append(xchlst)
            chlbls.append('Select the intensity data item')
            choices.append(ych)
            chlbls.append('Select the data item to be used for weighting')
            choices.append(such)
            chlbls.append('Divide intensities by data item')
            choices.append(['none']+modch)
            from .. import GSASIIctrlGUI as G2G
            res = G2G.MultipleChoicesSelector(choices,chlbls)
            if not res:
                self.errors = "Abort: data items not selected"
                return False
            xi,yi,sui,modi = res

            # now read in the values
            # x-values
            self.powderentry[0] = filename
            #self.powderentry[1] = pos # bank offset (N/A here)
            #self.powderentry[2] = 1 # xye file only has one bank
            self.idstring = os.path.basename(filename) + ': ' + blk
            if cf[blk].get('_diffrn_radiation_probe'):
                if cf[blk]['_diffrn_radiation_probe'] == 'neutron':
                    self.instdict['type'] = 'PNC'
                    #if cf[blk].get('_pd_meas_time_of_flight'): self.instdict['type'] = 'PNT' # not supported yet
                else:
                    self.instdict['type'] = 'PXC'
            if cf[blk].get('_diffrn_radiation_wavelength'):
                val = cf[blk]['_diffrn_radiation_wavelength']
                wl = []
                if type(val) is list:
                    for v in val:
                        w,e = cif.get_number_with_esd(v)
                        if w: wl.append(w)
                else:
                    w,e = cif.get_number_with_esd(val)
                    if w: wl.append(w)
                if wl:
                    if len(wl) > 1:
                        self.instdict['wave'] = wl
                    else:
                        self.instdict['wave'] = wl[0]
            if cf[blk].get('_diffrn_ambient_temperature'):
                val = cf[blk]['_diffrn_ambient_temperature']
                w,e = cif.get_number_with_esd(val)
                if w:
                    self.Sample['Temperature'] = w
        xcf = xch[xi]
        if type(xcf) is tuple:
            vals = [float(cf[blk][ixi]) for ixi in xcf]
            x = np.array([(i*vals[2] + vals[0]) for i in range(1 + int(0.5 + (vals[1]-vals[0])/vals[2]))])
        else:
            vl = []
            for val in cf[blk].get(xcf,'?'):
                v,e = cif.get_number_with_esd(val)
                if v is None: # not parsed
                    vl.append(np.nan)
                else:
                    vl.append(v)
            x = np.array(vl)
            if 'recip_len_q' in xcf and 'wave' in self.instdict:
                wl = self.instdict['wave']
                x = 2.*asind(wl*x/(4.*np.pi))
        # y-values
        ycf = ych[yi]
        vl = []
        v2 = []
        for val in cf[blk].get(ycf,'?'):
            v,e = cif.get_number_with_esd(val)
            if v is None: # not parsed
                vl.append(np.nan)
                v2.append(np.nan)
            else:
                vl.append(v)
                if e is None:
                    v2.append(np.sqrt(v))
                else:
                    v2.append(max(e,1.0))
        y = np.array(vl)
        w = 1./np.array(v2)**2
        # weights
        if sui == -1:
            # no weights
            vl = np.zeros(len(x)) + 1.
        else:
            vl = []
            sucf = such[sui]
            if sucf ==  '_pd_proc_ls_weight':
                for val in cf[blk].get(sucf,'?'):
                    v,e = cif.get_number_with_esd(val)
                    if v is None: # not parsed
                        vl.append(0.)
                    else:
                        vl.append(v)
            elif sucf ==  '_pd_proc_intensity_total':
                for val in cf[blk].get(sucf,'?'):
                    v,e = cif.get_number_with_esd(val)
                    if v is None: # not parsed
                        vl.append(0.)
                    elif v <= 0:
                        vl.append(1.)
                    else:
                        vl.append(1./v)
            elif sucf ==  '_pd_meas_counts_total':
                for val in cf[blk].get(sucf,'?'):
                    v,e = cif.get_number_with_esd(val)
                    if v is None: # not parsed
                        vl.append(0.)
                    elif v <= 0:
                        vl.append(1.)
                    else:
                        vl.append(1./v)
            else:
                for val in cf[blk].get(sucf,'?'):
                    v,e = cif.get_number_with_esd(val)
                    if v is None or e is None: # not parsed or no ESD
                        vl.append(np.nan)
                    elif e <= 0:
                        vl.append(1.)
                    else:
                        vl.append(1./(e*e))
#            w = np.array(vl)
        # intensity modification factor
        if modi >= 1:
            ycf = modch[modi-1]
            vl = []
            for val in cf[blk].get(ycf,'?'):
                v,e = cif.get_number_with_esd(val)
                if v is None: # not parsed
                    vl.append(np.nan)
                else:
                    vl.append(v)
            y /= np.array(vl)
            w /= np.array(vl)
        N = len(x)
        print ("CIF file, read from selected block")

        self.errors = "Error while storing read values"
        self.powderdata = [
                np.array(x), # x-axis values
                np.array(y), # powder pattern intensities
                np.array(w), # 1/sig(intensity)^2 values (weights)
                np.zeros(N), # calc. intensities (zero)
                np.zeros(N), # calc. background (zero)
                np.zeros(N), # obs-calc profiles
            ]
        return True
