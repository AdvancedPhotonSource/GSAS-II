# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2sfact_CIF: CIF import*
-----------------------------------
Read structure factors from a CIF reflection table.

'''
# routines to read in structure factors from a CIF
# 
import sys
import numpy as np
import os.path
import GSASIIIO as G2IO
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import CifFile as cif # PyCifRW from James Hester

class CIFhklReader(G2IO.ImportStructFactor):
    'Routines to import Phase information from a CIF file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist = ('.CIF','.cif','.FCF','.fcf','.HKL','.hkl'),
            strictExtension = False,
            formatName = 'CIF',
            longFormatName = 'CIF format structure factor file (.cif or .hkl)'
            )
    # Validate the contents
    def ContentsValidator(self, filepointer):
        'Use standard CIF validator'
        return self.CIFValidator(filepointer)

    def Reader(self, filename, filepointer, ParentFrame=None, **kwarg):
        '''Read single crystal data from a CIF.
        If multiple datasets are requested, use self.repeat and buffer caching.
        '''
        hklitems = [('_refln_index_h','_refln_index_k','_refln_index_l'),
                    ('_refln.index_h','_refln.index_k','_refln.index_l')]
        cellitems = [
            ('_cell_length_a','_cell_length_b','_cell_length_c',
             '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',),
            ('_cell.length_a','_cell.length_b','_cell.length_c',
             '_cell.angle_alpha','_cell.angle_beta','_cell.angle_gamma',),]

        Fdatanames = ('_refln_f_meas','_refln.f_meas','_refln.f_meas_au',
                      )
        
        F2datanames = ('_refln_f_squared_meas','_refln.f_squared_meas',
            '_refln_intensity_meas','_refln.intensity_meas',
                      )

        Idatanames = ('_refln_intensity_meas','_refln.intensity_meas',
                      ) # not used yet

        Isignames = ('_refln_intensity_meas_sigma','_refln.intensity_meas_sigma',
                      ) # not used yet

        Fcalcnames = ('_refln_f_calc','_refln.f_calc','_refln.f_calc_au',
                      )
        
        F2calcnames = ('_refln_f_squared_calc','_refln.f_squared_calc',
                      )

        Fsignames = ('_refln_f_meas_sigma','_refln.f_meas_sigma','_refln.f_meas_sigma_au',
                    '_refln_f_sigma',
                      )
        
        F2signames = ('_refln_f_squared_meas_sigma','_refln.f_squared_meas_sigma',
                      '_refln_f_squared_sigma',
                      '_refln_intensity_meas_sigma','_refln.intensity_meas_sigma',
                      '_refln.intensity_sigma',)

        phasenames = ('_refln_phase_calc','_refln.phase_calc',
                      )


        SGdataname = ('_symmetry_space_group_name_H-M', '_symmetry.space_group_name_H-M')
                     
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
                print 'Reusing previously parsed CIF'
            if cf is None:
                self.ShowBusy() # this can take a while
                try:
                    cf = G2IO.ReadCIF(filename)
                except Exception as detail:
                    self.errors = "Parse or reading of file failed in pyCifRW; check syntax of file in enCIFer or CheckCIF"
                    return False
                finally:
                    self.DoneBusy()
            # scan blocks for reflections
            self.errors = 'Error during scan of blocks for datasets'
            blklist = []
            for blk in cf.keys(): # scan for reflections, F or F2 values and cell lengths.
                # Ignore blocks that do not have structure factors and a cell
                blkkeys = [k.lower() for k in cf[blk].keys()]
                gotFo = False
                gotFo2 = False
                for i in range(2):
                    if hklitems[i][0] in blkkeys and hklitems[i][1] in blkkeys and hklitems[i][2] in blkkeys:
                        dnIndex = i
                        break
                else:
                    break # no reflections
                for dn in Fdatanames: 
                    if dn in blkkeys:
                        blklist.append(blk)
                        gotFo = True
                        break
                if gotFo: break
                for dn in F2datanames: 
                    if dn in blkkeys:
                        blklist.append(blk)
                        break
                else:
                    break
            if not blklist:
                selblk = None # no block to choose
            elif len(blklist) == 1: # only one choice
                selblk = 0
            else:                       # choose from options
                choice = []
                for blknm in blklist:
                    choice.append('')
                    # accumulate some info about this phase
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
                    for i,key in enumerate(cellitems[dnIndex]):
                        if i == 3: fmt = "%.f,"
                        if i == 5: fmt = "%.f"
                        val = cf[blknm].get(key)
                        if val is None: break
                        s += fmt % cif.get_number_with_esd(val)[0]
                    if s: choice[-1] += ', cell: ' + s
                    for dn in SGdataname:
                        sg = cf[blknm].get(dn)
                        if sg: 
                            choice[-1] += ', (' + sg.strip() + ')'
                            break
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
            self.errors = 'Error during reading of selected block'
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
            self.errors = 'Error during reading of reflections'
            # read in reflections
            try:
                refloop = blk.GetLoop(hklitems[0][0])
                dnIndex = 0
            except KeyError:
                try:
                    refloop = blk.GetLoop(hklitems[1][0])
                    dnIndex = 1
                except KeyError:
                    self.errors += "\nUnexpected: '_refln[-.]index_h not found!"
                    return False
            itemkeys = {}
            # prepare an index to the CIF reflection loop
            for i,key in enumerate(refloop.keys()):
                itemkeys[key.lower()] = i
                
            # scan for data names:
            F2dn = None
            Fdn = None
            F2cdn = None
            Fcdn = None
            F2sdn = None
            Fsdn = None
            Phdn = None
            FcalcPresent = False
            for dn in F2datanames:
                if dn in itemkeys:
                    F2dn = dn
                    for dn in F2calcnames:
                        if dn in itemkeys:
                            F2cdn = dn
                            FcalcPresent = True
                            break
                    for dn in F2signames:
                        if dn in itemkeys:
                            F2sdn = dn
                            break
                    break
            else:
                for dn in Fdatanames:
                    if dn in itemkeys:
                        Fdn = dn
                        for dn in Fcalcnames:
                            if dn in itemkeys:
                                Fcdn = dn
                                FcalcPresent = True
                                break
                        for dn in Fsignames:
                            if dn in itemkeys:
                                Fsdn = dn
                                break
                        break
                else:
                    msg = "\nno F or F2 loop value found in file\n"
                    msg += "A CIF reflection file needs to have at least one of\n"
                    for dn in F2datanames+Fdatanames:
                        msg += dn + ', '
                    self.errors += msg                        
                    return False
            for dn in phasenames:
                if dn in itemkeys:
                    Phdn = dn
                    break
                
            # loop over all reflections
            for item in refloop:
                F2c = 0.0
                sigF2 = 0.0
                HKL = []
                try:
                    for i in hklitems[dnIndex]: # '_refln[._]index_[hkl]'
                        num = itemkeys.get(i)
                        try:
                            HKL.append(int(item[num]))
                        except:
                            HKL.append('.')
                    #h,k,l,m,dsp,Fo2,sig,Fc2,Fot2,Fct2,phase,Ext
                    ref = HKL+[0,0,0,0,0, 0,0,0,0,0, 0] 
                    if F2dn:
                        F2 = item[itemkeys[F2dn]]
                        if '(' in F2:
                            F2, sigF2 = cif.get_number_with_esd(F2)
                            F2 = float(F2)
                            sigF2 = float(sigF2)
                        elif F2sdn:
                            F2 = float(F2)
                            sigF2 = float(item[itemkeys[F2sdn]])
                        else:
                            F2 = float(F2)
                        try:
                            if F2cdn:
                                F2c = float(item[itemkeys[F2cdn]])
                        except:
                            pass
                    else:
                        F = item[itemkeys[Fdn]]
                        if '(' in F:
                            F, sig = cif.get_number_with_esd(F)
                        elif Fsdn:
                            F = float(F)
                            sig = float(item[itemkeys[Fsdn]])
                        else:
                            F = float(F)
                            sig = 0.0
                        F2 = F**2
                        sigF2 = 2.0*F*sig
                        try:
                            if Fcdn:
                                Fc = float(item[itemkeys[Fcdn]])
                                F2c = Fc*Fc
                        except:
                            pass
                                
                    ref[8] = F2
                    ref[5] = F2
                    ref[6] = sigF2
                    ref[9] = F2c
                    ref[7] = F2c
                    try:
                        if Phdn:
                            ref[10] = float(item[itemkeys[Phdn]])
                    except:
                        pass
                except:
                    continue # skip over incompletely parsed reflections
                self.RefDict['RefList'].append(ref)
                self.RefDict['FF'].append({})
            self.RefDict['RefList'] = np.array(self.RefDict['RefList'])
            self.errors = 'Error during reading of dataset parameters'
            if blk.get('_diffrn_radiation_probe'):
                if blk['_diffrn_radiation_probe'] == 'neutron':
                    Type = 'SNC'
            elif blk.get('_diffrn_radiation.probe'):
                if blk['_diffrn_radiation.probe'] == 'neutron':
                    Type = 'SNC'
            else:
                type = 'SXC'
            self.RefDict['Type'] = Type
            if blk.get('_diffrn_radiation_wavelength'):
                wave = float(blk['_diffrn_radiation_wavelength'])
            elif blk.get('_diffrn_radiation.wavelength'):
                wave = float(blk['_diffrn_radiation.wavelength'])
            else:
                wave = 1.5418
            self.UpdateParameters(Type=Type,Wave=wave) # histogram type
            return True
        except Exception as detail:
            self.errors += '\n  '+str(detail)
            print self.formatName+' read error:'+str(detail) # for testing
            import traceback
            traceback.print_exc(file=sys.stdout)
        return False
