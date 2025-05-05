# -*- coding: utf-8 -*-
'''
'''
# Routines to import Phase information from CIF files
from __future__ import division, print_function
import sys
import random as ran
import numpy as np
import re
import copy
import os.path
from .. import GSASIIpath
from .. import GSASIIobj as G2obj
from .. import GSASIIspc as G2spc
from .. import GSASIIElem as G2elem
from .. import GSASIIlattice as G2lat
from .. import GSASIIfiles as G2fil
try:
    import CifFile as cif # PyCifRW from James Hester as a package
except ImportError:
    try:
        from .. import CifFile as cif # PyCifRW, as distributed w/G2 (old)
    except ImportError:
        cif = None
debug = GSASIIpath.GetConfigValue('debug')
#debug = False

class CIFPhaseReader(G2obj.ImportPhase):
    'Implements a phase importer from a possibly multi-block CIF file'
    def __init__(self):
        if cif is None:
            self.UseReader = False
            msg = 'CIFPhase Reader skipped because PyCifRW (CifFile) module is not installed.'
            G2fil.ImportErrorMsg(msg,{'CIF Phase importer':['pycifrw']})
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.CIF','.cif','.mcif'),
            strictExtension=False,
            formatName = 'CIF',
            longFormatName = 'Crystallographic Information File import'
            )

    def ContentsValidator(self, filename):
        fp = open(filename,'r')
        ok = self.CIFValidator(fp)
        #print('validator: ',ok)
        fp.close()
        return ok

    def Reader(self,filename, ParentFrame=None, usedRanIdList=[], **unused):
        if cif is None: # unexpected, but worth a specific error message
            print('Attempting to read a CIF without PyCifRW installed')
            raise Exception('Attempting to read a CIF without PyCifRW installed')
        isodistort_warnings = '' # errors that would prevent an isodistort analysis
        self.Phase = G2obj.SetNewPhase(Name='new phase') # create a new empty phase dict
        # make sure the ranId is really unique!
        while self.Phase['ranId'] in usedRanIdList:
            self.Phase['ranId'] = ran.randint(0,sys.maxsize)
        self.MPhase = G2obj.SetNewPhase(Name='new phase') # create a new empty phase dict
        # make sure the ranId is really unique!
        while self.MPhase['ranId'] in usedRanIdList:
            self.MPhase['ranId'] = ran.randint(0,sys.maxsize)
        returnstat = False
        cellitems = (
            '_cell_length_a','_cell_length_b','_cell_length_c',
            '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',)
#        cellwaveitems = (
#            '_cell_wave_vector_seq_id',
#            '_cell_wave_vector_x','_cell_wave_vector_y','_cell_wave_vector_z')
        reqitems = (
             '_atom_site_fract_x',
             '_atom_site_fract_y',
             '_atom_site_fract_z',
            )
        phasenamefields = (
            '_chemical_name_common',
            '_pd_phase_name',
            '_chemical_formula_sum'
            )
        try:
            cf = G2obj.ReadCIF(filename)
        except cif.StarError as msg:
            msg  = f'\nThis file does not have valid CIF syntax. Web site https://checkcif.iucr.org/ can help find CIF errors. If VESTA or https://addie.ornl.gov/conf_viewer is able to read this CIF, it may allow you to rewrite it as a valid file. \n\nError from PyCifRW: {msg}'
            self.errors = msg
            self.warnings += msg
            return False
        # scan blocks for structural info
        self.errors = 'Error during scan of blocks for datasets'
        str_blklist = []
        for blk in cf.keys():
            for r in reqitems+cellitems:
                if r not in cf[blk].keys():
                    break
            else:
                str_blklist.append(blk)
        if not str_blklist:
            selblk = None # no block to choose
        elif len(str_blklist) == 1: # only one choice
            selblk = 0
        else:                       # choose from options
            choice = []
            for blknm in str_blklist:
                choice.append('')
                # accumumlate some info about this phase
                choice[-1] += blknm + ': '
                for i in phasenamefields: # get a name for the phase
                    name = cf[blknm].get(i,'phase name').strip()
                    if name is None or name == '?' or name == '.':
                        continue
                    else:
                        choice[-1] += name.strip() + ', '
                        break
                na = len(cf[blknm].get("_atom_site_fract_x"))
                if na == 1:
                    choice[-1] += '1 atom'
                else:
                    choice[-1] += ('%d' % na) + ' atoms'
                choice[-1] += ', cell: '
                fmt = "%.2f,"
                for i,key in enumerate(cellitems):
                    if i == 3: fmt = "%.f,"
                    if i == 5: fmt = "%.f"
                    choice[-1] += fmt % cif.get_number_with_esd(
                        cf[blknm].get(key))[0]
                sg = cf[blknm].get("_symmetry_space_group_name_H-M",'')
                if not sg: sg = cf[blknm].get("_space_group_name_H-M_alt",'')
                if not sg: sg = cf[blknm].get("_space_group_ssg_name",'')
                if not sg: sg = cf[blknm].get("_space_group.magn_ssg_name_BNS",'')
                if not sg: sg = cf[blknm].get("_space_group.magn_ssg_name",'')
                #how about checking for super/magnetic ones as well? - reject 'X'?
                sg = sg.replace('_','')
                if sg: choice[-1] += ', (' + sg.strip() + ')'
            try:
                from .. import GSASIIctrlGUI as G2G
                selblk = G2G.PhaseSelector(choice,ParentFrame=ParentFrame,
                    title= 'Select a phase from one the CIF data_ blocks below',size=(600,100))
            except: # no wxPython
                selblk = 0
        self.errors = 'Error during reading of selected block'
#process selected phase
        if selblk is None:
            returnstat = False # no block selected or available
        else:   #do space group symbol & phase type first
            blknm = str_blklist[selblk]
            blk = cf[str_blklist[selblk]]
            E = True
            Super = False
            magnetic = False
            moddim = int(blk.get("_cell_modulation_dimension",'0'))
            if moddim:      #incommensurate
                if moddim > 1:
                    msg = 'more than 3+1 super symmetry is not allowed in GSAS-II'
                    self.errors = msg
                    self.warnings += '\n'+msg
                    return False
                if blk.get('_cell_subsystems_number'):
                    msg = 'Composite super structures not allowed in GSAS-II'
                    self.errors = msg
                    self.warnings += '\n'+msg
                    return False
                sspgrp = blk.get("_space_group_ssg_name",'')
                if not sspgrp:          #might be incommensurate magnetic
                    MSSpGrp = blk.get("_space_group.magn_ssg_name_BNS",'')
                    if not MSSpGrp:
                        MSSpGrp = blk.get("_space_group.magn_ssg_name",'')
                    if not MSSpGrp:
                        msg = 'No incommensurate space group name was found in the CIF.'
                        self.errors = msg
                        self.warnings += '\n'+msg
                        return False
                    if 'X' in MSSpGrp:
                        msg = 'Ad hoc incommensurate magnetic space group '+MSSpGrp+' is not allowed in GSAS-II'
                        self.warnings += '\n'+msg
                        self.errors = msg
                        return False
                    magnetic = True
                if 'X' in sspgrp:
                    msg = 'Ad hoc incommensurate space group '+sspgrp+' is not allowed in GSAS-II'
                    self.warnings += '\n'+msg
                    self.errors = msg
                    return False
                Super = True
                if magnetic:
                    sspgrp = MSSpGrp.split('(')
                    sspgrp[1] = "("+sspgrp[1]
                    SpGrp = G2spc.StandardizeSpcName(sspgrp[0])
                    if "1'" in SpGrp: sspgrp[1] = sspgrp[1][:-1]  #take off extra 's'; gets put back later
                    MSpGrp = sspgrp[0]
                    self.MPhase['General']['Type'] = 'magnetic'
                    self.MPhase['General']['AtomPtrs'] = [3,1,10,12]
                else:
                    sspgrp = sspgrp.split('(')
                    sspgrp[1] = "("+sspgrp[1]
                    SpGrp = sspgrp[0]
                    if "1'" in SpGrp:       #nonmagnetics can't be gray
                        SpGrp = SpGrp.replace("1'",'')
                        sspgrp[1] = sspgrp[1][:-1]  #take off extra 's'
#                    SpGrp = G2spc.StandardizeSpcName(SpGrp)
                    self.Phase['General']['Type'] = 'nuclear'
                if not SpGrp:
                    print (sspgrp)
                    self.warnings += 'Space group name '+sspgrp[0]+sspgrp[1]+' not recognized by GSAS-II'
                    return False
                SuperSg = sspgrp[1].replace('\\','').replace(',','')
                SuperVec = [[0,0,.1],False,4]
            else:   #not incommensurate
                SpGrp = blk.get("_symmetry_space_group_name_H-M",'')
                if not SpGrp:
                    SpGrp = blk.get("_space_group_name_H-M_alt",'')
                try:
                    SpGrp = G2spc.spgbyNum[int(blk.get('_symmetry_Int_Tables_number'))]
                except:
                    pass
                if not SpGrp:   #try magnetic
                    MSpGrp = blk.get("_space_group.magn_name_BNS",'')
                    if not MSpGrp:
                        MSpGrp = blk.get("_space_group_magn.name_BNS",'')
                        # if not MSpGrp:
                        #     msg = 'No recognizable space group name was found in the CIF.'
                        #     self.errors = msg
                        #     self.warnings += '\n'+msg
                        #     return False
                    SpGrp = blk.get('_parent_space_group.name_H-M_alt')
                    if not SpGrp:
                        SpGrp = blk.get('_parent_space_group.name_H-M')
                    if SpGrp and MSpGrp:
#                    SpGrp = MSpGrp.replace("'",'')
                        SpGrp = SpGrp[:2]+SpGrp[2:].replace('_','')   #get rid of screw '_'
                        if '_' in SpGrp[1]: SpGrp = SpGrp.split('_')[0]+SpGrp[3:]
                        SpGrp = G2spc.StandardizeSpcName(SpGrp)
                        magnetic = True
                        self.MPhase['General']['Type'] = 'magnetic'
                        self.MPhase['General']['AtomPtrs'] = [3,1,10,12]
                    elif not SpGrp:
                        print (MSpGrp)
                        self.warnings += 'No space group name was found in the CIF.'
                        #return False
                        SpGrp = 'P 1'
                else:
                    SpGrp = SpGrp.replace('_','').split('(')[0]
                    SpGrp = G2spc.fullHM2shortHM(SpGrp)
                    self.Phase['General']['Type'] = 'nuclear'
#process space group symbol
            E,SGData = G2spc.SpcGroup(SpGrp)
            if E and SpGrp:
                SpGrpNorm = G2spc.StandardizeSpcName(SpGrp)
                if SpGrpNorm:
                    E,SGData = G2spc.SpcGroup(SpGrpNorm)
            # if E:   #try lookup from number  - found full symbol?
            #     SpGrpNorm = G2spc.spgbyNum[int(blk.get('_symmetry_Int_Tables_number'))]
            #     if SpGrpNorm:
            #         E,SGData = G2spc.SpcGroup(SpGrpNorm)
            # nope, try the space group "out of the Box"
            if E:
                self.warnings += 'ERROR in space group symbol '+SpGrp
                self.warnings += '\nThe space group has been set to "P 1". '
                self.warnings += "Change this in phase's General tab."
                self.warnings += '\nAre there spaces separating axial fields?\n\nError msg: '
                self.warnings += G2spc.SGErrors(E)
                SGData = G2obj.P1SGData # P 1
            self.Phase['General']['SGData'] = SGData
            # save symmetry operators, if specified (use to check for origin 1 in GSASIIdataGUI OnImportPhase)
            ops = blk.get("_symmetry_equiv_pos_as_xyz")   # try older name 1st
            if ops:
                self.SymOps['xyz'] = ops
            elif blk.get("_space_group_symop_operation_xyz"):
                self.SymOps['xyz'] = blk.get("_space_group_symop_operation_xyz")
            else:
                if 'xyz' in self.SymOps: del self.SymOps['xyz']
            if Super:
                E,SSGData = G2spc.SSpcGroup(SGData,SuperSg)
                if E:
                    self.warnings += 'Invalid super symmetry symbol '+SpGrp+SuperSg
                    self.warnings += '\n'+E
                    SuperSg = SuperSg[:SuperSg.index(')')+1]
                    self.warnings += '\nNew super symmetry symbol '+SpGrp+SuperSg
                    E,SSGData = G2spc.SSpcGroup(SGData,SuperSg)
                self.Phase['General']['SSGData'] = SSGData
                if magnetic:
                    self.MPhase['General']['SGData'] = SGData
                    self.MPhase['General']['SGData']['MagSpGrp'] = MSSpGrp.replace(',','').replace('\\','')
                    self.MPhase['General']['SSGData'] = SSGData

            if magnetic:    #replace std operaors with those from cif file - probably not the same!
                SGData['SGFixed'] = True
                SGData['SGOps'] = []
                SGData['SGCen'] = []
                if Super:
                    SSGData['SSGOps'] = []
                    SSGData['SSGCen'] = []
                    try:
                        sgoploop = blk.GetLoop('_space_group_symop_magn_ssg_operation.id')
                        sgcenloop = blk.GetLoop('_space_group_symop_magn_ssg_centering.id')
                        opid = sgoploop.GetItemPosition('_space_group_symop_magn_ssg_operation.algebraic')[1]
                        centid = sgcenloop.GetItemPosition('_space_group_symop_magn_ssg_centering.algebraic')[1]
                    except KeyError:        #old mag cif names
                        sgoploop = blk.GetLoop('_space_group_symop.magn_ssg_id')
                        sgcenloop = blk.GetLoop('_space_group_symop.magn_ssg_centering_id')
                        opid = sgoploop.GetItemPosition('_space_group_symop.magn_ssg_operation_algebraic')[1]
                        centid = sgcenloop.GetItemPosition('_space_group_symop.magn_ssg_centering_algebraic')[1]
                    spnflp = []
                    for op in sgoploop:
                        M,T,S = G2spc.MagSSText2MTS(op[opid])
                        SSGData['SSGOps'].append([np.array(M,dtype=float),T])
                        SGData['SGOps'].append([np.array(M,dtype=float)[:3,:3],T[:3]])
                        spnflp.append(S)
                    censpn = []
                    for cent in sgcenloop:
                        M,C,S = G2spc.MagSSText2MTS(cent[centid])
                        SSGData['SSGCen'].append(C)
                        SGData['SGCen'].append(C[:3])
                        censpn += list(np.array(spnflp)*S)
                    self.MPhase['General']['SSGData'] = SSGData
                else:
                    try:
                        sgoploop = blk.GetLoop('_space_group_symop_magn_operation.id')
                        opid = sgoploop.GetItemPosition('_space_group_symop_magn_operation.xyz')[1]
                        try:
                            sgcenloop = blk.GetLoop('_space_group_symop_magn_centering.id')
                            centid = sgcenloop.GetItemPosition('_space_group_symop_magn_centering.xyz')[1]
                        except KeyError:
                            sgcenloop = None
                    except KeyError:
                        try:
                            sgoploop = blk.GetLoop('_space_group_symop_magn.id')
                            sgcenloop = blk.GetLoop('_space_group_symop_magn_centering.id')
                            opid = sgoploop.GetItemPosition('_space_group_symop_magn_operation.xyz')[1]
                            centid = sgcenloop.GetItemPosition('_space_group_symop_magn_centering.xyz')[1]
                        except KeyError:        #old mag cif names
                            sgoploop = blk.GetLoop('_space_group_symop.magn_id')
                            sgcenloop = blk.GetLoop('_space_group_symop.magn_centering_id')
                            opid = sgoploop.GetItemPosition('_space_group_symop.magn_operation_xyz')[1]
                            centid = sgcenloop.GetItemPosition('_space_group_symop.magn_centering_xyz')[1]
                    spnflp = []
                    for op in sgoploop:
                        try:
                            M,T,S = G2spc.MagText2MTS(op[opid])
                            SGData['SGOps'].append([np.array(M,dtype=float),T])
                            spnflp.append(S)
                        except KeyError:
                            self.warnings += 'Space group operator '+op[opid]+' is not recognized by GSAS-II'
                            return False
                    censpn = []
                    if sgcenloop:
                        for cent in sgcenloop:
                            M,C,S = G2spc.MagText2MTS(cent[centid])
                            SGData['SGCen'].append(C)
                            censpn += list(np.array(spnflp)*S)
                    else:
                            M,C,S = G2spc.MagText2MTS('x,y,z,+1')
                            SGData['SGCen'].append(C)
                            censpn += list(np.array(spnflp)*S)
                self.MPhase['General']['SGData'] = SGData
                self.MPhase['General']['SGData']['SpnFlp'] = censpn
                G2spc.GenMagOps(SGData)         #set magMom
                self.MPhase['General']['SGData']['MagSpGrp'] = MSpGrp
                if "1'" in MSpGrp:
                    SGData['SGGray'] = True
                MagPtGp = blk.get('_space_group.magn_point_group')
                if not MagPtGp:
                    MagPtGp = blk.get('_space_group.magn_point_group_name')
                if not MagPtGp:
                    MagPtGp = blk.get('_space_group_magn.point_group_name')
                self.MPhase['General']['SGData']['MagPtGp'] = MagPtGp

            # cell parameters
            cell = []
            for lbl in cellitems:
                cell.append(cif.get_number_with_esd(blk[lbl])[0])
            Volume = G2lat.calc_V(G2lat.cell2A(cell))
            self.Phase['General']['Cell'] = [False,]+cell+[Volume,]
            if magnetic:
                self.MPhase['General']['Cell'] = [False,]+cell+[Volume,]
            if Super:
                waveloop = blk.GetLoop('_cell_wave_vector_seq_id')
                waveDict = dict(waveloop.items())
                SuperVec = [[cif.get_number_with_esd(waveDict['_cell_wave_vector_x'][0].replace('?','0'))[0],
                    cif.get_number_with_esd(waveDict['_cell_wave_vector_y'][0].replace('?','0'))[0],
                    cif.get_number_with_esd(waveDict['_cell_wave_vector_z'][0].replace('?','0'))[0]],False,4]

            # read in atoms
            self.errors = 'Error during reading of atoms'
            atomlbllist = [] # table to look up atom IDs
            atomloop = blk.GetLoop('_atom_site_label')
            atomkeys = [i.lower() for i in atomloop.keys()]
            if not blk.get('_atom_site_type_symbol'):
                isodistort_warnings += '\natom types are missing. \n Check & revise atom types as needed'
            if magnetic:
                try:
                    magmoment = '_atom_site_moment.label'
                    magatomloop = blk.GetLoop(magmoment)
                    magatomkeys = [i.lower() for i in magatomloop.keys()]
                    magatomlabels = blk.get(magmoment)
                    G2MagDict = {'_atom_site_moment.label': 0,
                                 '_atom_site_moment.crystalaxis_x':7,
                                 '_atom_site_moment.crystalaxis_y':8,
                                 '_atom_site_moment.crystalaxis_z':9}
                except KeyError:
                    magmoment = '_atom_site_moment_label'
                    magatomloop = blk.GetLoop(magmoment)
                    magatomkeys = [i.lower() for i in magatomloop.keys()]
                    magatomlabels = blk.get(magmoment)
                    G2MagDict = {'_atom_site_moment_label': 0,
                                 '_atom_site_moment_crystalaxis_x':7,
                                 '_atom_site_moment_crystalaxis_y':8,
                                 '_atom_site_moment_crystalaxis_z':9}

            if blk.get('_atom_site_aniso_label'):
                anisoloop = blk.GetLoop('_atom_site_aniso_label')
                anisokeys = [i.lower() for i in anisoloop.keys()]
                anisolabels = blk.get('_atom_site_aniso_label')
            else:
                anisoloop = None
                anisokeys = []
                anisolabels = []
            if Super:
#                occFloop = None
                occCloop = None
#                occFdict = {}
                occCdict = {}
#                displSloop = None
                displFloop = None
                MagFloop = None
#                displSdict = {}
                displFdict = {}
                MagFdict = {}
                UijFloop = None
                UijFdict = {}
                #occupancy modulation
#                if blk.get('_atom_site_occ_Fourier_atom_site_label'):
#                    occFloop = blk.GetLoop('_atom_site_occ_Fourier_atom_site_label')
#                    occFdict = dict(occFloop.items())
                if blk.get('_atom_site_occ_special_func_atom_site_label'):  #Crenel (i.e. Block Wave) occ
                    occCloop = blk.GetLoop('_atom_site_occ_special_func_atom_site_label')
                    occCdict = dict(occCloop.items())
                #position modulation
                if blk.get('_atom_site_displace_Fourier_atom_site_label'):
                    displFloop = blk.GetLoop('_atom_site_displace_Fourier_atom_site_label')
                    displFdict = dict(displFloop.items())
#                if blk.get('_atom_site_displace_special_func_atom_site_label'): #sawtooth
#                    displSloop = blk.GetLoop('_atom_site_displace_special_func_atom_site_label')
#                    displSdict = dict(displSloop.items())
                #U modulation
                if blk.get('_atom_site_U_Fourier_atom_site_label'):
                    UijFloop = blk.GetLoop('_atom_site_U_Fourier_atom_site_label')
                    UijFdict = dict(UijFloop.items())
                #Mag moment modulation
                if blk.get('_atom_site_moment_Fourier_atom_site_label'):
                    MagFloop = blk.GetLoop('_atom_site_moment_Fourier_atom_site_label')
                    MagFdict = dict(MagFloop.items())
                    Mnames =  ['_atom_site_moment_fourier_atom_site_label',
                               '_atom_site_moment_fourier_axis','_atom_site_moment_fourier_wave_vector_seq_id',
                               '_atom_site_moment_fourier_param_sin','_atom_site_moment_fourier_param_cos']
                elif blk.get('_atom_site_moment_Fourier.atom_site_label'):
                    MagFloop = blk.GetLoop('_atom_site_moment_Fourier.atom_site_label')
                    MagFdict = dict(MagFloop.items())
                    Mnames =  ['_atom_site_moment_fourier.atom_site_label',
                               '_atom_site_moment_fourier.axis','_atom_site_moment_fourier.wave_vector_seq_id',
                               '_atom_site_moment_fourier_param.sin','_atom_site_moment_fourier_param.cos']
            self.Phase['Atoms'] = []
            if magnetic:
                self.MPhase['Atoms'] = []
            G2AtomDict = {  '_atom_site_type_symbol' : 1,
                            '_atom_site_label' : 0,
                            '_atom_site_fract_x' : 3,
                            '_atom_site_fract_y' : 4,
                            '_atom_site_fract_z' : 5,
                            '_atom_site_occupancy' : 6,
                            '_atom_site_aniso_u_11' : 11,
                            '_atom_site_aniso_u_22' : 12,
                            '_atom_site_aniso_u_33' : 13,
                            '_atom_site_aniso_u_12' : 14,
                            '_atom_site_aniso_u_13' : 15,
                            '_atom_site_aniso_u_23' : 16, }

            ranIdlookup = {}
            for aitem in atomloop:
                atomlist = ['','','',0.,0.,0.,1.0,'',0.,'I',0.01, 0.,0.,0.,0.,0.,0.,]
                for val,key in zip(aitem,atomkeys):
                    col = G2AtomDict.get(key,-1)
                    if col >= 3:
                        atomlist[col] = cif.get_number_with_esd(val)[0]
                        if col >= 11: atomlist[9] = 'A' # if any Aniso term is defined, set flag
                    elif col is not None and col != -1:
                        atomlist[col] = val
                    elif key in ('_atom_site_thermal_displace_type',
                               '_atom_site_adp_type'):   #Iso or Aniso?
                        if val.lower() == 'uani':
                            atomlist[9] = 'A'
                    elif key == '_atom_site_u_iso_or_equiv':
                        uisoval = cif.get_number_with_esd(val)[0]
                        if uisoval is not None:
                            atomlist[10] = uisoval
                    elif key == '_atom_site_b_iso_or_equiv':
                        uisoval = cif.get_number_with_esd(val)[0]
                        if uisoval is not None:
                            atomlist[10] = uisoval/(8*np.pi**2)
                if not atomlist[1] and atomlist[0]:
                    typ = atomlist[0].rstrip('0123456789-+')
                    if G2elem.CheckElement(typ):
                        atomlist[1] = typ
                    if not atomlist[1]:
                        atomlist[1] = 'Xe'
                        self.warnings += ' Atom type '+typ+' not recognized; Xe assumed\n'
                if atomlist[0] in anisolabels: # does this atom have aniso values in separate loop?
                    atomlist[9] = 'A'
                    for val,key in zip(anisoloop.GetKeyedPacket('_atom_site_aniso_label',atomlist[0]),anisokeys):
                        col = G2AtomDict.get(key)
                        if col:
                            atomlist[col] = cif.get_number_with_esd(val)[0]
                if None in atomlist[11:17]:
                    atomlist[9] = 'I'
                    atomlist[11:17] =  [0.,0.,0.,0.,0.,0.]
                atomlist[7],atomlist[8] = G2spc.SytSym(atomlist[3:6],SGData)[:2]
                atomlist[1] = G2elem.FixValence(atomlist[1])
                atomlist.append(ran.randint(0,sys.maxsize)) # add a random Id
                self.Phase['Atoms'].append(atomlist)
                ranIdlookup[atomlist[0]] = atomlist[-1]
                if atomlist[0] in atomlbllist:
                    self.warnings += ' ERROR: repeated atom label: '+atomlist[0]
                else:
                    atomlbllist.append(atomlist[0])

                if magnetic and atomlist[0] in magatomlabels:
                    matomlist = atomlist[:7]+[0.,0.,0.,]+atomlist[7:]
                    for mval,mkey in zip(magatomloop.GetKeyedPacket(magmoment,atomlist[0]),magatomkeys):
                        mcol = G2MagDict.get(mkey,-1)
                        if mcol > 0:
                            matomlist[mcol] = cif.get_number_with_esd(mval)[0]
                    self.MPhase['Atoms'].append(matomlist)
                if Super:
                    Sfrac = np.zeros((4,2))
                    Sadp = np.zeros((4,12))
                    Spos = np.zeros((4,6))
                    Smag = np.zeros((4,6))
                    nim = -1
                    waveType = 'Fourier'
                    if  occCdict:
                        for i,item in enumerate(occCdict['_atom_site_occ_special_func_atom_site_label']):
                            if item == atomlist[0]:
                                waveType = 'Crenel'
                                val = occCdict['_atom_site_occ_special_func_crenel_c'][i]
                                Sfrac[0][0] = cif.get_number_with_esd(val)[0]
                                val = occCdict['_atom_site_occ_special_func_crenel_w'][i]
                                Sfrac[0][1] = cif.get_number_with_esd(val)[0]
                                nim = 1

                    if nim >= 0:
                        Sfrac = [waveType,]+[[sfrac,False] for sfrac in Sfrac[:nim]]
                    else:
                        Sfrac = []
                    nim = -1
                    if displFdict:
                        for i,item in enumerate(displFdict['_atom_site_displace_fourier_atom_site_label']):
                            if item == atomlist[0]:
                                waveType = 'Fourier'
                                ix = ['x','y','z'].index(displFdict['_atom_site_displace_fourier_axis'][i])
                                im = int(displFdict['_atom_site_displace_fourier_wave_vector_seq_id'][i])
                                if im != nim:
                                    nim = im
                                val = displFdict['_atom_site_displace_fourier_param_sin'][i]
                                Spos[im-1][ix] = cif.get_number_with_esd(val)[0]
                                val = displFdict['_atom_site_displace_fourier_param_cos'][i]
                                Spos[im-1][ix+3] = cif.get_number_with_esd(val)[0]
                    if nim >= 0:
                        Spos = [waveType,]+[[spos,False] for spos in Spos[:nim]]
                    else:
                        Spos = []
                    nim = -1
                    if UijFdict:
                        nim = -1
                        for i,item in enumerate(UijFdict['_atom_site_u_fourier_atom_site_label']):
                            if item == atomlist[0]:
                                ix = ['U11','U22','U33','U12','U13','U23'].index(UijFdict['_atom_site_u_fourier_tens_elem'][i])
                                im = int(UijFdict['_atom_site_u_fourier_wave_vector_seq_id'][i])
                                if im != nim:
                                    nim = im
                                val = UijFdict['_atom_site_u_fourier_param_sin'][i]
                                Sadp[im-1][ix] = cif.get_number_with_esd(val)[0]
                                val = UijFdict['_atom_site_u_fourier_param_cos'][i]
                                Sadp[im-1][ix+6] = cif.get_number_with_esd(val)[0]
                    if nim >= 0:
                        Sadp = ['Fourier',]+[[sadp,False] for sadp in Sadp[:nim]]
                    else:
                        Sadp = []
                    nim = -1
                    if MagFdict:
                        nim = -1
                        for i,item in enumerate(MagFdict[Mnames[0]]):
                            if item == atomlist[0]:
                                ix = ['x','y','z'].index(MagFdict[Mnames[1]][i])
                                im = int(MagFdict[Mnames[2]][i])
                                if im != nim:
                                    nim = im
                                val = MagFdict[Mnames[3]][i]
                                Smag[im-1][ix] = cif.get_number_with_esd(val)[0]
                                val = MagFdict[Mnames[4]][i]
                                Smag[im-1][ix+3] = cif.get_number_with_esd(val)[0]
                    if nim >= 0:
                        Smag = ['Fourier',]+[[smag,False] for smag in Smag[:nim]]
                    else:
                        Smag = []
                    SSdict = {'SS1':{'Sfrac':Sfrac,'Spos':Spos,'Sadp':Sadp,'Smag':Smag}}
                    if magnetic and atomlist[0] in magatomlabels:
                        matomlist.append(SSdict)
                    atomlist.append(SSdict)
            if len(atomlbllist) != len(self.Phase['Atoms']):
                isodistort_warnings += '\nRepeated atom labels prevents ISODISTORT decode'
            for lbl in phasenamefields: # get a name for the phase
                try:
                    name = blk.get(lbl)
                except:
                    continue
                if name is None:
                    continue
                name = name.strip()
                if name == '?' or name == '.':
                    continue
                else:
                    break
            else: # no name found, use block name for lack of a better choice; for isodistort use filename
                name = blknm
            if 'isodistort' in name:
                self.Phase['General']['Name'] = os.path.split(filename)[1].split('.')[0]
            else:
                self.Phase['General']['Name'] = name.strip()
            self.Phase['General']['Super'] = Super
            self.Phase = copy.deepcopy(self.Phase)      #clean copy
            if magnetic:
                self.MPhase = copy.deepcopy(self.MPhase)    #clean copy
                self.MPhase['General']['Type'] = 'magnetic'
                self.MPhase['General']['Name'] = name.strip()+' mag'
                self.MPhase['General']['Super'] = Super
                if Super:
                    self.MPhase['General']['Modulated'] = True
                    self.MPhase['General']['SuperVec'] = SuperVec
                    self.MPhase['General']['SuperSg'] = SuperSg
                    if self.MPhase['General']['SGData']['SGGray']:
                        self.MPhase['General']['SuperSg'] += 's'
                if 'mcif' not in filename:
                    self.Phase = copy.deepcopy(self.MPhase)
                    del self.MPhase
            else:
                self.MPhase = None
            if Super:
                if self.Phase['General']['SGData']['SGGray']:
                    SGData = self.Phase['General']['SGData']
                    SGData['SGGray'] = False
                    ncen = len(SGData['SGCen'])
                    SGData['SGCen'] = SGData['SGCen'][:ncen//2]
                    self.Phase['General']['SGData'].update(SGData)
                self.Phase['General']['Modulated'] = True
                self.Phase['General']['SuperVec'] = SuperVec
                self.Phase['General']['SuperSg'] = SuperSg
                if not self.Phase['General']['SGData']['SGFixed']:
                    self.Phase['General']['SSGData'] = G2spc.SSpcGroup(SGData,SuperSg)[1]
            if self.ISODISTORT_test(blk):
                if isodistort_warnings:
                    self.warnings += isodistort_warnings
                else:
                    self.errors = "Error while processing ISODISTORT constraints"
                    self.ISODISTORT_proc(blk,atomlbllist,ranIdlookup,filename)
                    self.errors = ""
            returnstat = True
        if self.SymOps.get('xyz',[]): # did the phase supply symmetry operations?
            # check if they agree with GSAS-II's setting
            good = G2spc.CompareSym(self.SymOps['xyz'],
                                     SGData=self.Phase['General']['SGData'])
            if not good:
                centro = False
                if '-x,-y,-z' in [i.replace(' ','').lower() for i in self.SymOps['xyz']]:
                    centro = True
                msg = 'This CIF uses a space group setting not compatible with GSAS-II.\n'
                if self.Phase['General']['SGData']['SpGrp'] in G2spc.spg2origins:
                    msg += '\nThis is likely due to the space group being set in Origin 1 rather than 2.\n'
                msg += '''
Do you want to use Bilbao's "CIF to Standard Setting" web service to
transform this into a standard setting?
'''
                if self.Phase['General']['SGData']['SpGrp'] in G2spc.spg2origins and not centro:
                    msg += '''
If you say "no" here, later, you will get the chance to apply a simple origin
shift later as an alternative to the above.'''
                else:
                    msg += '\nIf you say "no" here you will need to do this yourself manually.'
                try:
                    from .. import GSASIIctrlGUI as G2G
                    ans = G2G.askQuestion(ParentFrame,msg,'xform structure?')
                except Exception as err: # fails if non-interactive (no wxPython)
                    print(err)
                    print('\nCIF symops do not agree with GSAS-II, calling Bilbao "CIF to Standard Setting" web service.\n')
                    ans = True
                if ans:
                    from .. import SUBGROUPS
                    SUBGROUPS.createStdSetting(filename,self)
        return returnstat

    def ISODISTORT_test(self,blk):
        '''Test if there is any ISODISTORT information in CIF

        At present only _iso_displacivemode... and _iso_occupancymode... are
        tested.
        '''
        for i in ('_iso_displacivemode_label',
                  '_iso_occupancymode_label'):
            if blk.get(i): return True
        return False

    def ISODISTORT_proc(self,blk,atomlbllist,ranIdlookup,filename):
        '''Process ISODISTORT items to create constraints etc.
        Constraints are generated from information extracted from
        loops beginning with _iso_ and are placed into
        self.Constraints, which contains a list of
        :ref:`constraints tree items <Constraint_definitions_table>`
        and one dict.
        The dict contains help text for each generated ISODISTORT variable

        At present only _iso_displacivemode... and _iso_occupancymode... are
        processed. Not yet processed: _iso_magneticmode...,
        _iso_rotationalmode... & _iso_strainmode...
        '''
        varLookup = {'dx':'dAx','dy':'dAy','dz':'dAz','do':'Afrac'}
        'Maps ISODISTORT parm names to GSAS-II names'
        # used for all types of modes
        self.Constraints = []
        explaination = {}
        G2obj.AddPhase2Index(self,filename)   # put phase info into Var index
        #----------------------------------------------------------------------
        # read in the ISODISTORT displacement modes
        #----------------------------------------------------------------------
        if blk.get('_iso_displacivemode_label'):
            modelist = []
            shortmodelist = []
            modedispl = []
            idlist = []
            for id,lbl,val in zip(
                blk.get('_iso_displacivemode_ID'),
                blk.get('_iso_displacivemode_label'),
                blk.get('_iso_displacivemode_value')):
                idlist.append(int(id))
                modelist.append(lbl)
                modedispl.append(float(val))
                ISODISTORT_shortLbl(lbl,shortmodelist) # shorten & make unique
            # just in case the items are not ordered increasing by id, sort them here
            modelist = [i for i,j in sorted(zip(modelist,idlist),key=lambda k:k[1])]
            shortmodelist = [i for i,j in sorted(zip(shortmodelist,idlist),key=lambda k:k[1])]
            # read in the coordinate offset variables names and map them to G2 names/objects
            coordVarLbl = []
            G2varObj = []
            G2coordOffset = []
            idlist = []
            error = False
            ParentCoordinates = {}
            for lbl,exp in zip(
                blk.get('_iso_coordinate_label'),
                blk.get('_iso_coordinate_formula') ):
                if '_' in lbl:
                    albl = lbl[:lbl.rfind('_')]
                    vlbl = lbl[lbl.rfind('_')+1:]
                else:
                    self.warnings += ' ERROR: _iso_coordinate_label not parsed: '+lbl
                    error = True
                    continue
                if vlbl not in 'xyz' or len(vlbl) != 1:
                    self.warnings += ' ERROR: _iso_coordinate_label coordinate not parsed: '+lbl
                    error = True
                    continue
                i = 'xyz'.index(vlbl)
                if not ParentCoordinates.get(albl):
                    ParentCoordinates[albl] = [None,None,None]
                if '+' in exp:
                    val = exp.split('+')[0].strip()
                    val = G2fil.FormulaEval(val)
                elif '-' in exp:
                    val = exp.split('-')[0].strip()
                    val = G2fil.FormulaEval(val)
                else:
                    val = G2fil.FormulaEval(exp)
                if val is None:
                    self.warnings += ' ERROR: _iso_coordinate_formula coordinate not interpreted: '+lbl
                    error = True
                    continue
                else:
                    ParentCoordinates[albl][i] = val
            if error:
                print (self.warnings)
                raise Exception("Error decoding variable labels")

            error = False
            coordOffset = {xyz:i for i,xyz in enumerate(('dx','dy','dz'))}
            for id,lbl,val in zip(
                blk.get('_iso_deltacoordinate_ID'),
                blk.get('_iso_deltacoordinate_label'),
                blk.get('_iso_deltacoordinate_value') ):
                idlist.append(int(id))
                coordVarLbl.append(lbl)
                if '_' in lbl:
                    albl = lbl[:lbl.rfind('_')]
                    vlbl = lbl[lbl.rfind('_')+1:]
                else:
                    self.warnings += ' ERROR: _iso_deltacoordinate_label not parsed: '+lbl
                    error = True
                    continue
                if albl not in atomlbllist:
                    self.warnings += ' ERROR: _iso_deltacoordinate_label atom not found: '+lbl
                    error = True
                    continue
                # else:
                #     anum = atomlbllist.index(albl)
                var = varLookup.get(vlbl)
                if not var:
                    self.warnings += ' ERROR: _iso_deltacoordinate_label variable not found: '+lbl
                    error = True
                    continue
                G2coordOffset.append(ParentCoordinates[albl][coordOffset[vlbl]] - float(val))
                G2varObj.append(G2obj.G2VarObj(
                    (self.Phase['ranId'],None,var,ranIdlookup[albl])
                    ))
            if error:
                raise Exception("Error decoding variable labels")
            # just in case the items are not ordered increasing by id, sort them here
            coordVarLbl = [i for i,j in sorted(zip(coordVarLbl,idlist),key=lambda k:k[1])]
            G2varObj = [i for i,j in sorted(zip(G2varObj,idlist),key=lambda k:k[1])]

            if len(G2varObj) != len(modelist):
                print ("non-square input")
                raise Exception("Rank of _iso_displacivemode != _iso_deltacoordinate")

            # normalization constants
            normlist = []
            idlist = []
            for id,exp in zip(
                blk.get('_iso_displacivemodenorm_ID'),
                blk.get('_iso_displacivemodenorm_value'),
                ):
                idlist.append(int(id))
                normlist.append(float(exp))
            normlist = [i for i,j in sorted(zip(normlist,idlist),key=lambda k:k[1])]
            # get mapping of modes to atomic coordinate displacements
            displacivemodematrix = np.zeros((len(G2varObj),len(G2varObj)))
            for row,col,val in zip(
                blk.get('_iso_displacivemodematrix_row'),
                blk.get('_iso_displacivemodematrix_col'),
                blk.get('_iso_displacivemodematrix_value'),):
                displacivemodematrix[int(row)-1,int(col)-1] = float(val)
            # Invert to get mapping of atom displacements to modes
            Var2ModeMatrix = np.linalg.inv(displacivemodematrix)
            # create the constraints
            modeVarList = []
            for i,(row,norm) in enumerate(zip(Var2ModeMatrix,normlist)):
                constraint = []
                for j,(lbl,k) in enumerate(zip(coordVarLbl,row)):
                    if k == 0: continue
                    constraint.append([k/norm,G2varObj[j]])
                modeVar = G2obj.G2VarObj(
                    (self.Phase['ranId'],None,shortmodelist[i],None))
                modeVarList.append(modeVar)
                constraint += [modeVar,False,'f']
                self.Constraints.append(constraint)
            #----------------------------------------------------------------------
            # save the ISODISTORT info for "mode analysis"
            if 'ISODISTORT' not in self.Phase: self.Phase['ISODISTORT'] = {}
            self.Phase['ISODISTORT'].update({
                # coordinate items
                'IsoVarList' : coordVarLbl,
                'G2VarList' : G2varObj,
                'G2coordOffset' : G2coordOffset,
                'G2parentCoords' : [ParentCoordinates[item] for item in ParentCoordinates], #Assumes python 3.7 dict ordering!
                # mode items
                'IsoModeList' : modelist,
                'G2ModeList' : modeVarList,
                'NormList' : normlist,
                'modeDispl' : modedispl,
                'ISOmodeDispl' : copy.deepcopy(modedispl),
                # transform matrices
                'Var2ModeMatrix' : Var2ModeMatrix,
                'Mode2VarMatrix' : displacivemodematrix,
                })
            # make explaination dictionary
            for mode,shortmode in zip(modelist,shortmodelist):
                modeVar = G2obj.G2VarObj(
                    (self.Phase['ranId'],None,shortmode,None))
                explaination[modeVar] = ("Full ISODISTORT name for " +
                                            shortmode + " is " + str(mode))
        #----------------------------------------------------------------------
        # now read in the ISODISTORT occupancy modes
        #----------------------------------------------------------------------
        if blk.get('_iso_occupancymode_label'):
            modelist = []
            shortmodelist = []
            idlist = []
            for id,lbl in zip(
                blk.get('_iso_occupancymode_ID'),
                blk.get('_iso_occupancymode_label')):
                idlist.append(int(id))
                modelist.append(lbl)
                ISODISTORT_shortLbl(lbl,shortmodelist) # shorten & make unique
            # just in case the items are not ordered increasing by id, sort them here
            modelist = [i for i,j in sorted(zip(modelist,idlist),key=lambda k:k[1])]
            shortmodelist = [i for i,j in sorted(zip(shortmodelist,idlist),key=lambda k:k[1])]
            # read in the coordinate offset variables names and map them to G2 names/objects
            occVarLbl = []
            G2varObj = []
            idlist = []
            error = False
            for id,lbl in zip(
                blk.get('_iso_deltaoccupancy_ID'),
                blk.get('_iso_deltaoccupancy_label') ):
                idlist.append(int(id))
                occVarLbl.append(lbl)
                if '_' in lbl:
                    albl = lbl[:lbl.rfind('_')]
                    vlbl = lbl[lbl.rfind('_')+1:]
                else:
                    self.warnings += ' ERROR: _iso_deltaoccupancy_label not parsed: '+lbl
                    error = True
                    continue
                if albl not in atomlbllist:
                    self.warnings += ' ERROR: _iso_deltaoccupancy_label atom not found: '+lbl
                    error = True
                    continue
                # else:
                #     anum = atomlbllist.index(albl)
                var = varLookup.get(vlbl)
                if not var:
                    self.warnings += ' ERROR: _iso_deltaoccupancy_label variable not found: '+lbl
                    error = True
                    continue
                G2varObj.append(G2obj.G2VarObj(
                    (self.Phase['ranId'],None,var,ranIdlookup[albl])
                    ))
            if error:
                raise Exception("Error decoding variable labels")
            # just in case the items are not ordered increasing by id, sort them here
            occVarLbl = [i for i,j in sorted(zip(occVarLbl,idlist),key=lambda k:k[1])]
            G2varObj = [i for i,j in sorted(zip(G2varObj,idlist),key=lambda k:k[1])]

            if len(G2varObj) != len(modelist):
                print ("non-square input")
                raise Exception("Rank of _iso_occupancymode != _iso_deltaoccupancy")

            error = False
            ParentOcc = {}
            for lbl,exp in zip(
                blk.get('_iso_occupancy_label'),
                blk.get('_iso_occupancy_formula') ):
                if '_' in lbl:
                    albl = lbl[:lbl.rfind('_')]
                    vlbl = lbl[lbl.rfind('_')+1:]
                else:
                    self.warnings += ' ERROR: _iso_occupancy_label not parsed: '+lbl
                    error = True
                    continue
                if vlbl != 'occ':
                    self.warnings += ' ERROR: _iso_occupancy_label coordinate not parsed: '+lbl
                    error = True
                    continue
                if '+' in exp:
                    val = exp.split('+')[0].strip()
                    val = G2fil.FormulaEval(val)
                    if val is None:
                        self.warnings += ' ERROR: _iso_occupancy_formula coordinate not interpreted: '+lbl
                        error = True
                        continue
                    ParentOcc[albl] = val
            if error:
                raise Exception("Error decoding occupancy labels")
            # get mapping of modes to atomic coordinate displacements
            occupancymodematrix = np.zeros((len(G2varObj),len(G2varObj)))
            for row,col,val in zip(
                blk.get('_iso_occupancymodematrix_row'),
                blk.get('_iso_occupancymodematrix_col'),
                blk.get('_iso_occupancymodematrix_value'),):
                occupancymodematrix[int(row)-1,int(col)-1] = float(val)
            # Invert to get mapping of atom displacements to modes
            occupancymodeInvmatrix = np.linalg.inv(occupancymodematrix)
            # create the constraints
            modeVarList = []
            for i,row in enumerate(occupancymodeInvmatrix):
                constraint = []
                for j,(lbl,k) in enumerate(zip(occVarLbl,row)):
                    if k == 0: continue
                    constraint.append([k,G2varObj[j]])
                modeVar = G2obj.G2VarObj(
                    (self.Phase['ranId'],None,shortmodelist[i],None))
                modeVarList.append(modeVar)
                constraint += [modeVar,False,'f']
                self.Constraints.append(constraint)
            # normilization constants
            normlist = []
            idlist = []
            for id,exp in zip(
                blk.get('_iso_occupancymodenorm_ID'),
                blk.get('_iso_occupancymodenorm_value'),
                ):
                idlist.append(int(id))
                normlist.append(float(exp))
            normlist = [i for i,j in sorted(zip(normlist,idlist),key=lambda k:k[1])]
            #----------------------------------------------------------------------
            # save the ISODISTORT info for "mode analysis"
            if 'ISODISTORT' not in self.Phase: self.Phase['ISODISTORT'] = {}
            self.Phase['ISODISTORT'].update({
                # coordinate items
                'OccVarList' : occVarLbl,
                'G2OccVarList' : G2varObj,
                'BaseOcc' : ParentOcc,
                # mode items
                'OccModeList' : modelist,
                'G2OccModeList' : modeVarList,
                'OccNormList' : normlist,
                # transform matrices
                'Var2OccMatrix' : occupancymodeInvmatrix,
                'Occ2VarMatrix' : occupancymodematrix,
                })
            # make explaination dictionary
            for mode,shortmode in zip(modelist,shortmodelist):
                modeVar = G2obj.G2VarObj(
                    (self.Phase['ranId'],None,shortmode,None))
                explaination[modeVar] = ("Full ISODISTORT name for " +
                                            shortmode + " is " + str(mode))
        if explaination: self.Constraints.append(explaination)
        #----------------------------------------------------------------------
        # done with read
        #----------------------------------------------------------------------

        def fmtEqn(i,head,l,var,k):
            'format a section of a row of variables and multipliers'
            if np.isclose(k,0): return head,l
            if len(head) + len(l) > 65:
                print(head+l)
                head = 20*' '
                l = ''
            if k < 0 and i > 0:
                l += ' - '
                k = -k
            elif i > 0:
                l += ' + '
            if k == 1:
                l += '%s ' % str(var)
            else:
                l += '%.3f * %s' % (k,str(var))
            return head,l

        # debug: show displacive mode var to mode relations
        if debug and 'IsoVarList' in self.Phase['ISODISTORT']:
            print('\n' + 70*'=')
            print('ISO modes from Iso coordinate vars (using Var2ModeMatrix, IsoVarList, G2VarList & G2ModeList)' )
            for i,row in enumerate(self.Phase['ISODISTORT']['Var2ModeMatrix']):
                norm = self.Phase['ISODISTORT']['NormList'][i]
                head = '  ' + str(self.Phase['ISODISTORT']['G2ModeList'][i]) + ' = ('
                line = ''
                for j,(lbl,k) in enumerate(zip(coordVarLbl,row)):
                    var = self.Phase['ISODISTORT']['IsoVarList'][j]
                    head,line = fmtEqn(j,head,line,var,k)
                print(head+line+') / {:.3g}'.format(norm))
                head = '              = '
                line = ''
                for j,(lbl,k) in enumerate(zip(coordVarLbl,row)):
                    var = self.Phase['ISODISTORT']['IsoVarList'][j]
                    head,line = fmtEqn(j,head,line,var,k/norm)
                print(head+line)
            print('\nConstraints')
            for c in self.Constraints:
                if type(c) is dict: continue
                if c[-1] != 'f': continue
                line = ''
                head = '  ' + str(c[-3]) + ' = '
                for j,(k,var) in enumerate(c[:-3]):
                    head,line = fmtEqn(j,head,line,var,k)
                print(head+line)

            # Get the ISODISTORT offset values
            coordVarDelta = {}
            for lbl,val in zip(
                blk.get('_iso_deltacoordinate_label'),
                blk.get('_iso_deltacoordinate_value'),):
                coordVarDelta[lbl] = float(val)
            modeVarDelta = {}
            for lbl,val in zip(
                blk.get('_iso_displacivemode_label'),
                blk.get('_iso_displacivemode_value'),):
                modeVarDelta[lbl] = cif.get_number_with_esd(val)[0]

            print('\n' + 70*'=')
            print('Confirming inverse mode relations computed from displacement values',
                      '\nusing Var2ModeMatrix, NormList, IsoVarList')
            # compute the mode values from the reported coordinate deltas
            for i,(row,n) in enumerate(zip(self.Phase['ISODISTORT']['Var2ModeMatrix'],
                                           self.Phase['ISODISTORT']['NormList'])):
                line = ''
                print(str(self.Phase['ISODISTORT']['IsoModeList'][i])+' = ')
                head = '  = ('
                for j,(lbl,k) in enumerate(zip(coordVarLbl,row)):
                    head,line = fmtEqn(j,head,line,lbl,k)
                print(head+line+') / '+('%.3f'%n))
                line = ''
                head = '  = ('
                vsum = 0.
                for j,(lbl,k) in enumerate(zip(coordVarLbl,row)):
                    val = "{:3g}".format(coordVarDelta[lbl])
                    head,line = fmtEqn(j,head,line,val,k)
                    vsum += coordVarDelta[lbl] * k
                print(head+line+') / '+('%.3f'%n))
                fileval = modeVarDelta[self.Phase['ISODISTORT']['IsoModeList'][i]]
                print("{} = {:4g} (value read from CIF = {:4g})\n".format(
                    self.Phase['ISODISTORT']['IsoModeList'][i], vsum, fileval))

            print( 70*'=')
            print('Direct displacement relations computed from ISO modes in CIF',
                      '\nusing Mode2VarMatrix, NormList, IsoModeList, IsoVarList',)
            # compute the coordinate displacements from the reported mode values
            for lbl,row in zip(self.Phase['ISODISTORT']['IsoVarList'],
                               self.Phase['ISODISTORT']['Mode2VarMatrix']):
                l = ''
                s = 0.0
                head = lbl+' ='
                for j,(k,n) in enumerate(zip(row,self.Phase['ISODISTORT']['NormList'])):
                    if k == 0: continue
                    if len(l) > 65:
                        print(head,l)
                        head = 20*' '
                        l = ''
                    l1 = ''
                    k1 = k
                    if j > 0 and k < 0:
                        k1 = -k
                        l1 = ' - '
                    elif j > 0:
                        l1 += ' + '
                    l += '{:} {:3g} * {:4g} * {:}'.format(
                        l1, k1, n, self.Phase['ISODISTORT']['IsoModeList'][j])

                    s += n * modeVarDelta[self.Phase['ISODISTORT']['IsoModeList'][j]] * k
                print(head,l)
                print(lbl,'=',s)
                print(lbl,'==>',str(self.Phase['ISODISTORT']['G2VarList'][i]),'\n')
            DeltaCoords = {}
            for i,lbl,row in zip(range(len(coordVarLbl)),coordVarLbl,displacivemodematrix):
                s = 0.0
                for j,(k,n) in enumerate(zip(row,self.Phase['ISODISTORT']['NormList'])):
                    s += k * n * modeVarDelta[self.Phase['ISODISTORT']['IsoModeList'][j]]
                at,d = lbl.rsplit('_',1)
                if at not in DeltaCoords:
                    DeltaCoords[at] = [0,0,0]
                if d == 'dx':
                    DeltaCoords[at][0] = s
                elif d == 'dy':
                    DeltaCoords[at][1] = s
                elif d == 'dz':
                    DeltaCoords[at][2] = s
                #else:
                #    print('unexpected',d)

            print( 70*'=')
            print('Coordinate checks')
            print("\nxyz's Computed from ISO mode values, as above")
            for at in sorted(DeltaCoords):
                s = at
                for i in range(3):
                    s += '  '
                    s += str(ParentCoordinates[at][i]+DeltaCoords[at][i])
                print(s)

            # determine the coordinate delta values from deviations from the parent structure
            print("\nxyz Values read directly from CIF")
            for atmline in self.Phase['Atoms']:
                lbl = atmline[0]
                x,y,z = atmline[3:6]
                print( lbl,x,y,z)

            print('\n' + 70*'=')
            print("G2 short name ==> ISODISTORT full name",
                      " (from IsoModeList and G2ModeList)")
            for mode,G2mode in zip(self.Phase['ISODISTORT']['IsoModeList'],
                                   self.Phase['ISODISTORT']['G2ModeList']):
                print('{} ==> {}'.format(str(G2mode), mode))
            print('\nConstraint help dict info')
            for i in self.Constraints:
                if type(i) is dict:
                    for key in i:
                        print('\t',key,':',i[key])
            print( 70*'=')

        #======================================================================
        # debug: show occupancy mode var to mode relations
        if debug and 'G2OccVarList' in self.Phase['ISODISTORT']:
            # coordinate items
            #occVarLbl = self.Phase['ISODISTORT']['OccVarList']
            G2varObj = self.Phase['ISODISTORT']['G2OccVarList']
            #ParentOcc = self.Phase['ISODISTORT']['BaseOcc']
            # mode items
            modelist = self.Phase['ISODISTORT']['OccModeList']
            modeVarList = self.Phase['ISODISTORT']['G2OccModeList']
            normlist = self.Phase['ISODISTORT']['OccNormList']
            # transform matrices
            #occupancymodeInvmatrix = self.Phase['ISODISTORT']['Var2OccMatrix']
            #occupancymodematrix = self.Phase['ISODISTORT']['Occ2VarMatrix']

            print( 70*'=')
            print('\nVar2OccMatrix' ,'OccVarList' )
            for i,row in enumerate(occupancymodeInvmatrix):
                l = ''
                for j,(lbl,k) in enumerate(zip(occVarLbl,row)):
                    if k == 0: continue
                    if l: l += ' + '
                    #l += lbl+' * '+str(k)
                    l += str(G2varObj[j])+' * '+str(k)
                print( str(i) + ': '+str(modeVarList[i])+' = '+l)

            # Get the ISODISTORT offset values
            occVarDelta = {}
            for lbl,val in zip(
                blk.get('_iso_deltaoccupancy_label'),
                blk.get('_iso_deltaoccupancy_value'),):
                occVarDelta[lbl] = float(val)
            modeVarDelta = {}
            for lbl,val in zip(
                blk.get('_iso_occupancymode_label'),
                blk.get('_iso_occupancymode_value'),):
                modeVarDelta[lbl] = cif.get_number_with_esd(val)[0]

            print( 70*'=')
            print('\nInverse relations using Var2OccModeMatrix, OccNormList, OccVarList')
            # compute the mode values from the reported coordinate deltas
            for i,(row,n) in enumerate(zip(occupancymodeInvmatrix,normlist)):
                l = ''
                for lbl,k in zip(occVarLbl,row):
                    if k == 0: continue
                    if l: l += ' + '
                    l += lbl+' * '+str(k)
                print('a'+str(i)+' = '+str(modeVarList[i])+' = ('+l+')/'+str(n))
            print('\nCalculation checks\n')
            for i,(row,n) in enumerate(zip(occupancymodeInvmatrix,normlist)):
                #l = ''
                sl = ''
                s = 0.
                for lbl,k in zip(occVarLbl,row):
                    if k == 0: continue
                    #if l: l += ' + '
                    #l += lbl+' * '+str(k)
                    if sl: sl += ' + '
                    sl += str(occVarDelta[lbl])+' * '+str(k)
                    s += occVarDelta[lbl] * k
                print(str(modeVarList[i]),'=','('+sl+') / ',n,'=',s/n)
                print(' ?= ',modeVarDelta[modelist[i]])
                print()

            print( 70*'=')
            print('\nDirect relations using Occ2VarMatrix, OccNormList, OccVarList')
            # compute the coordinate displacements from the reported mode values
            Occ = {}
            for i,lbl,row in zip(range(len(occVarLbl)),occVarLbl,occupancymodematrix):
                l = ''
                s = 0.0
                for j,(k,n) in enumerate(zip(row,normlist)):
                    if k == 0: continue
                    if l: l += ' + '
                    l += str(n)+' * '+str(modeVarList[j])+' * '+str(k)
                    s += n * modeVarDelta[modelist[j]] * k
                print( lbl,'=',str(G2varObj[i]),'=',l,'=',s,'\n')
                j = lbl.split('_')[0]
                Occ[j] = ParentOcc[j]+s

            # determine the coordinate delta values from deviations from the parent structure
            print('\nOccupancy from CIF vs computed')
            for atmline in self.Phase['Atoms']:
                lbl = atmline[0]
                if lbl in Occ: print( lbl,atmline[6],Occ[lbl])

            print( 70*'=')
            print('\nGenerated constraints')
            for i in self.Constraints:
                if type(i) is dict:
                    print('\nconstraint help dict')
                    for key in i:
                        print('\t',key,':',i[key])
                elif i[-1] == 'f':
                    print('\n\t',i[-3],' =')
                    for m,j in i[:-3]:
                        print('\t\t+',m,' * ',j,'   ',repr(j))
                else:
                    print('  unexpected: ',repr(i))
            print("\nG2name ==> ISODISTORT full name",
                      ".Phase['ISODISTORT']['OccModeList']",
                      ".Phase['ISODISTORT']['G2OccModeList']")
            for mode,G2mode in zip(modelist,modeVarList):
                print("  ?::"+str(G2mode),' ==>', mode)

def ISODISTORT_shortLbl(lbl,shortmodelist):
    '''Shorten model labels and remove special characters
    '''
    regexp = re.match(r'.*?\[.*?\](.*?)\(.*?\)\[.*?\](.*)',lbl)
    # assume lbl is of form SSSSS[x,y,z]AAAA(a,b,...)[...]BBBBB
    # where SSSSS is the parent spacegroup, [x,y,z] is a location, etc.
    if regexp:
        # this extracts the AAAAA and BBBBB parts of the string
        lbl = regexp.expand(r'\1\2') # parse succeeded, make a short version
    else:
        # try form SSSSS[x,y,z]AAAA(a,b,...)BBBBB
        regexp = re.match(r'.*?\[.*?\](.*?)\(.*?\)(.*)',lbl)
        if regexp:
            lbl = regexp.expand(r'\1\2') # parse succeeded, make a short version
    lbl = lbl.replace('order','o')
    lbl = lbl.replace('+','_')
    lbl = lbl.replace('-','_')
    G2obj.MakeUniqueLabel(lbl,shortmodelist) # make unique and add to list
