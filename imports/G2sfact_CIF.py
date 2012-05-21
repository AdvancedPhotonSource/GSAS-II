# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2012-02-13 11:33:35 -0600 (Mon, 13 Feb 2012) $
# $Author: vondreele & toby $
# $Revision: 482 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/G2importphase.py $
# $Id: G2importphase.py 482 2012-02-13 17:33:35Z vondreele $
########### SVN repository information ###################
# routines to read in structure factors from a CIF
# 
import sys
import numpy as np
import os.path
import GSASIIIO as G2IO
import CifFile as cif # PyCifRW from James Hester
import urllib

class CIFhklReader(G2IO.ImportStructFactor):
    'Routines to import Phase information from a CIF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.CIF','.cif','.HKL','.hkl'),
            strictExtension=False,
            formatName = 'CIF',
            longFormatName = 'CIF format structure factor file (.cif or .hkl)'
            )
    # Validate the contents
    def ContentsValidator(self, filepointer):
        for i,line in enumerate(filepointer):
            if i >= 1000: break
            ''' Encountered only blank lines or comments in first 1000
            lines. This is unlikely, but assume it is CIF since we are
            even less likely to find a file with nothing but hashes and
            blank lines'''
            line = line.strip()
            if len(line) == 0:
                continue # ignore blank lines
            elif line.startswith('#'):
                continue # ignore comments
            elif line.startswith('data_'):
                return True
            else:
                return False # found something else
        return True
    def Reader(self,filename,filepointer, ParentFrame=None, **kwarg):
        hklitems = ('_refln_index_h','_refln_index_k','_refln_index_l')
        cellitems = (
            '_cell_length_a','_cell_length_b','_cell_length_c',
            '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',)
        phasenamefields = (
            '_chemical_name_common',
            '_pd_phase_name',
            '_chemical_formula_sum'
            )
        rdbuffer = kwarg.get('buffer')
        cf = None
        try:
            if self.repeat and rdbuffer is not None:
                cf = rdbuffer.get('lastcif')
                print 'Reuse previously parsed CIF'
            if cf is None:
                self.ShowBusy() # this can take a while
                ciffile = 'file:'+urllib.pathname2url(filename)
                cf = cif.ReadCif(ciffile)
                self.DoneBusy()
            # scan blocks for reflections
            blklist = []
            blktype = []
            for blk in cf.keys():
                blkkeys = [k.lower() for k in cf[blk].keys()]
                gotFo = True
                gotFo2 = True
                for r in hklitems:
                    if r not in blkkeys:
                        gotFo = False
                        gotFo2 = False
                if '_refln_f_squared_meas' not in blkkeys:
                    gotFo = False
                if '_refln_f_squared_meas' not in blkkeys:
                    gotFo2 = False
                if gotFo or gotFo2:
                    blklist.append(blk)
                    blktype.append(gotFo2)
            if not blklist:
                selblk = None # no block to choose
            elif len(blklist) == 1: # only one choice
                selblk = 0
            else:                       # choose from options
                choice = []
                for blknm in blklist:
                    choice.append('')
                    # accumumlate some info about this phase
                    choice[-1] += blknm + ': '
                    for i in phasenamefields: # get a name for the phase
                        name = cf[blknm].get(i)
                        if name is None or name == '?' or name == '.':
                            continue
                        else:
                            choice[-1] += name.strip()[:20] + ', '
                            break
                    s = ''
                    fmt = "%.2f,"
                    for i,key in enumerate(cellitems):
                        if i == 3: fmt = "%.f,"
                        if i == 5: fmt = "%.f"
                        val = cf[blknm].get(key)
                        if val is None: break
                        s += fmt % cif.get_number_with_esd(val)[0]
                    if s: choice[-1] += ', cell: ' + s
                    sg = cf[blknm].get("_symmetry_space_group_name_H-M")
                    if sg: choice[-1] += ', (' + sg.strip() + ')'
                choice.append('Import all of the above')
                if self.repeat: # we were called to repeat the read
                    selblk = self.repeatcount
                    self.repeatcount += 1
                    if self.repeatcount >= len(blklist): self.repeat = False
                else:
                    selblk = self.BlockSelector(
                        choice,
                        ParentFrame=ParentFrame,
                        title='Select a dataset from one the CIF data_ blocks below',
                        size=(600,100),
                        header='Dataset Selector')
            if selblk is None:
                return False # no block selected or available
            if selblk >= len(blklist): # all blocks selected
                selblk = 0
                self.repeat = True
                if rdbuffer is not None:
                    rdbuffer['lastcif'] = cf # save the parsed cif for the next loop
                self.repeatcount = 1
            blknm = blklist[selblk]
            blk = cf[blklist[selblk]]
            self.objname = os.path.basename(filename)+':'+str(blknm)
            # read in reflections
            refloop = blk.GetLoop('_refln_index_h')
            itemkeys = {}
            # prepare an index to the CIF reflection loop
            for i,key in enumerate(refloop.keys()):
                itemkeys[key.lower()] = i
            if '_refln_f_squared_calc' in itemkeys:
                FcalcPresent = True
            elif '_refln_f_calc' in itemkeys:
                FcalcPresent = True
            else:
                FcalcPresent = False
            for item in refloop:
                HKL = []
                for i in hklitems: # ('_refln_index_h','_refln_index_k','_refln_index_l')
                    num = itemkeys.get(i)
                    try:
                        HKL.append(int(item[num]))
                    except:
                        HKL.append('.')
                #h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,...
                ref = HKL+[0,0,0,0,0,0,0,0,[],[],0,{}] 
                if '_refln_f_squared_meas' in itemkeys:
                    try:
                        Fsq = float(item[itemkeys['_refln_f_squared_meas']])
                        ref[8] = Fsq
                        ref[5] = Fsq
                    except:
                        try:
                            Fsq,Esd = item[itemkeys['_refln_f_squared_meas']].split('(')
                            Fsq = float(Fsq)
                            ref[8] = Fsq
                            ref[5] = Fsq
                            ref[6] = float(Esd[:-1])
                        except:
                            pass
                    if  '_refln_f_squared_sigma' in itemkeys:
                        try:
                            ref[6] = float(item[itemkeys['_refln_f_squared_sigma']])
                        except:
                            pass                            
                elif '_refln_f_meas' in itemkeys:
                    try:
                        Fsq = float(item[itemkeys['_refln_f_meas']])**2
                        ref[8] = Fsq
                        ref[5] = Fsq
                    except:
                        try:
                            Fsq,Esd = item[itemkeys['_refln_f_squared_meas']].split('(')
                            Fsq = float(Fsq)**2
                            ref[8] = Fsq
                            ref[5] = Fsq
                            ref[6] = 2.0*sqrt(Fsq)*float(Esd[:-1])
                        except:
                            pass                                
                    if  '_refln_f_sigma' in itemkeys:
                        try:
                            ref[6] = 2.*sqrt(ref[8])*float(item[itemkeys['_refln_f_sigma']])
                        except:
                            pass                                
                if '_refln_f_squared_calc' in itemkeys:
                    try:
                        Fsq = float(item[itemkeys['_refln_f_squared_calc']])
                        ref[9] = Fsq
                        ref[7] = Fsq
                    except:
                        pass                                
                elif '_refln_f_calc' in itemkeys:
                    try:
                        Fsq = float(item[itemkeys['_refln_f_calc']])**2
                        ref[9] = Fsq
                        ref[7] = Fsq
                    except:
                        pass                                
                if '_refln_phase_calc' in itemkeys:
                    try:
                        ref[10] = float(item[itemkeys['_refln_phase_calc']])
                    except:
                        pass                                
                self.RefList.append(ref)
            self.UpdateControls(Type='Fosq',FcalcPresent=FcalcPresent) # set Fobs type & if Fcalc values are loaded
            if blk.get('_diffrn_radiation_probe'):
                if blk['_diffrn_radiation_probe'] == 'neutron':
                    type = 'SNC'
            else:
                type = 'SXC'
            if blk.get('_diffrn_radiation_wavelength'):
                wave = blk['_diffrn_radiation_wavelength']
            else:
                wave = None
            self.UpdateParameters(Type=type,Wave=wave) # histogram type
            return True
        except Exception as detail:
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
        return False
