########### SVN repository information ###################
# $Date: 2012-02-13 11:33:35 -0600 (Mon, 13 Feb 2012) $
# $Author: vondreele & toby $
# $Revision: 482 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/G2importphase_CIF.py $
# $Id: G2importphase_CIF.py 482 2012-02-13 17:33:35Z vondreele $
########### SVN repository information ###################
# Routines to import Phase information from CIF files
import sys
import random as ran
import urllib
import numpy as np
import GSASIIIO as G2IO
import GSASIIobj as G2obj
import GSASIIspc as G2spc
import GSASIIlattice as G2lat
import GSASIIpy3 as G2p3
import CifFile as cif # PyCifRW from James Hester

class ISODISTORTPhaseReader(G2IO.ImportPhase):
    varLookup = {'dx':'dAx','dy':'dAy','dz':'dAz','do':'Afrac'}
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.CIF','.cif','.txt'),
            strictExtension=False,
            formatName = 'ISODISTORT CIF',
            longFormatName = 'CIF from ISODISTORT import'
            )
    def ContentsValidator(self, filepointer):
        filepointer.seek(0) # rewind the file pointer
        lines = filepointer.read(10000) # scan the first 10K bytes
        # quick tests for an ISODISTORT output file
        if lines.find('ISODISTORT') != -1:
            return True
        if lines.find('data_isodistort-') != -1:
            return True
        if lines.find('_iso_') != -1:
            return True
        return False
    def Reader(self,filename,filepointer, ParentFrame=None, usedRanIdList=[], **unused):
        import re
        self.Phase = G2IO.SetNewPhase(Name='new phase',SGData=G2IO.P1SGData) # create a new empty phase dict
        # make sure the ranId is really unique!
        while self.Phase['ranId'] in usedRanIdList:
            self.Phase['ranId'] = ran.randint(0,sys.maxint)
        returnstat = False
        cellitems = (
            '_cell_length_a','_cell_length_b','_cell_length_c',
            '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',)
        reqitems = (
             '_atom_site_type_symbol',
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
            self.ShowBusy() # this can take a while
            ciffile = 'file:'+urllib.pathname2url(filename)
            cf = cif.ReadCif(ciffile)
            # scan blocks for structural info
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
                        name = cf[blknm].get(i).strip()
                        if name is None or name == '?' or name == '.':
                            continue
                        else:
                            choice[-1] += name.strip()[:20] + ', '
                            break
                    na = len(cf[blknm].get("_atom_site_fract_x"))
                    if na == 1:
                        choice[-1] += '1 atom'
                    else:
                        choice[-1] += ('%d' % nd) + ' atoms'
                    choice[-1] += ', cell: '
                    fmt = "%.2f,"
                    for i,key in enumerate(cellitems):
                        if i == 3: fmt = "%.f,"
                        if i == 5: fmt = "%.f"
                        choice[-1] += fmt % cif.get_number_with_esd(
                            cf[blknm].get(key))[0]
                    sg = cf[blknm].get("_symmetry_space_group_name_H-M")
                    if sg: choice[-1] += ', (' + sg.strip() + ')'
                selblk = self.PhaseSelector(
                    choice,
                    ParentFrame=ParentFrame,
                    title= 'Select a phase from one the CIF data_ blocks below',
                    size=(600,100)
                    )
            if selblk is None:
                returnstat = False # no block selected or available
            else:
                blknm = str_blklist[selblk]
                blk = cf[str_blklist[selblk]]
                SpGrp = blk.get("_symmetry_space_group_name_H-M")
                if SpGrp:
                     E,SGData = G2spc.SpcGroup(SpGrp)
                if E:
                    self.warnings += ' ERROR in space group symbol '+SpGrp
                    self.warnings += ' N.B.: make sure spaces separate axial fields in symbol' 
                    self.warnings += G2spc.SGErrors(E)
                else:
                    self.Phase['General']['SGData'] = SGData
                # cell parameters
                cell = []
                for lbl in cellitems:
                    cell.append(cif.get_number_with_esd(blk[lbl])[0])
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
                self.Phase['General']['Cell'] = [False,]+cell+[Volume,]
                # read in atoms
                atomlbllist = [] # table to look up atom IDs
                atomloop = blk.GetLoop('_atom_site_label')
                atomkeys = [i.lower() for i in atomloop.keys()]
                if blk.get('_atom_site_aniso_label'):
                    anisoloop = blk.GetLoop('_atom_site_aniso_label')
                    anisokeys = [i.lower() for i in anisoloop.keys()]
                else:
                    anisoloop = None
                    anisokeys = []
                self.Phase['Atoms'] = []
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
                    atomlist = ['','','',0,0,0,1.0,'',0,'I',0.01,0,0,0,0,0,0,0]
                    atomlist[-1] = ran.randint(0,sys.maxint) # add a random Id
                    while atomlist[-1] in ranIdlookup:
                        atomlist[-1] = ran.randint(0,sys.maxint) # make it unique                        
                    for val,key in zip(aitem,atomkeys):
                        col = G2AtomDict.get(key)
                        if col >= 3:
                            atomlist[col] = cif.get_number_with_esd(val)[0]
                        elif col is not None:
                            atomlist[col] = val
                        elif key in ('_atom_site_thermal_displace_type',
                                   '_atom_site_adp_type'):   #Iso or Aniso?
                            if val.lower() == 'uani':
                                atomlist[9] = 'A'
                        elif key == '_atom_site_u_iso_or_equiv':
                            atomlist[10] =cif.get_number_with_esd(val)[0]
                    ulbl = '_atom_site_aniso_label'
                    if  atomlist[9] == 'A' and atomlist[0] in blk.get(ulbl):
                        for val,key in zip(anisoloop.GetKeyedPacket(ulbl,atomlist[0]),
                                           anisokeys):
                            col = G2AtomDict.get(key)
                            if col:
                                atomlist[col] = cif.get_number_with_esd(val)[0]
                    atomlist[7],atomlist[8] = G2spc.SytSym(atomlist[3:6],SGData)
                    self.Phase['Atoms'].append(atomlist)
                    ranIdlookup[atomlist[0]] = atomlist[-1]
                    if atomlist[0] in atomlbllist:
                        self.warnings += ' ERROR: repeated atom label: '+atomlist[0]
                        returnstat = False # cannot use this
                    else:
                        atomlbllist.append(atomlist[0])
                if len(atomlbllist) != len(self.Phase['Atoms']):
                    raise Exception,'Repeated atom labels prevents ISODISTORT decode'
                for lbl in phasenamefields: # get a name for the phase
                    name = blk.get(lbl)
                    if name is None:
                        continue
                    name = name.strip()
                    if name == '?' or name == '.':
                        continue
                    else:
                        break
                else: # no name found, use block name for lack of a better choice
                    name = blknm
                self.Phase['General']['Name'] = name.strip()[:20]
                #----------------------------------------------------------------------
                # now read in the ISODISTORT displacement modes
                #----------------------------------------------------------------------
                self.Constraints = []
                explaination = {}
                if blk.get('_iso_displacivemode_label'):
                    modelist = []
                    shortmodelist = []
                    for lbl in blk.get('_iso_displacivemode_label'):
                        modelist.append(lbl)
                        # assume lbl is of form SSSSS[x,y,z]AAAA(a,b,...)BBBBB
                        # where SSSSS is the parent spacegroup, [x,y,z] is a location 
                        regexp = re.match(r'.*?\[.*?\](.*?)\(.*?\)(.*)',lbl)
                        # this extracts the AAAAA and BBBBB parts of the string
                        if regexp:
                            lbl = regexp.expand(r'\1\2') # parse succeeded, make a short version
                        G2obj.MakeUniqueLabel(lbl,shortmodelist) # make unique and add to list
                    # read in the coordinate offset variables names and map them to G2 names/objects
                    coordVarLbl = []
                    G2varLbl = []
                    G2varObj = []
                    error = False
                    for lbl in blk.get('_iso_deltacoordinate_label'):
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
                        else:
                            anum = atomlbllist.index(albl)
                        var = self.varLookup.get(vlbl)
                        if not var:
                            self.warnings += ' ERROR: _iso_deltacoordinate_label variable not found: '+lbl
                            error = True
                            continue
                        G2varLbl.append('::'+var+':'+str(anum)) # variable name, less phase ID
                        G2varObj.append(G2obj.G2VarObj(
                            (self.Phase['ranId'],None,var,ranIdlookup[albl])
                            ))
                    if error:
                        raise Exception,"Error decoding variable labels"

                    if len(G2varObj) != len(modelist):
                        print "non-square input"
                        raise Exception,"Rank of _iso_displacivemode != _iso_deltacoordinate"

                    error = False
                    ParentCoordinates = {}
                    for lbl,exp in zip(
                        blk.get('_iso_coordinate_label'),
                        blk.get('_iso_coordinate_formula'),
                        ):
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
                            val = G2p3.FormulaEval(val)
                            if val is None:
                                self.warnings += ' ERROR: _iso_coordinate_formula coordinate not interpreted: '+lbl
                                error = True
                                continue
                        ParentCoordinates[albl][i] = val
                    if error:
                        print self.warnings
                        raise Exception,"Error decoding variable labels"
                    # get mapping of modes to atomic coordinate displacements
                    displacivemodematrix = np.zeros((len(G2varObj),len(G2varObj)))
                    for row,col,val in zip(
                        blk.get('_iso_displacivemodematrix_row'),
                        blk.get('_iso_displacivemodematrix_col'),
                        blk.get('_iso_displacivemodematrix_value'),):
                        displacivemodematrix[int(row)-1,int(col)-1] = float(val)
                    # Invert to get mapping of atom displacements to modes
                    displacivemodeInvmatrix = np.linalg.inv(displacivemodematrix)
                    # create the constraints
                    for i,row in enumerate(displacivemodeInvmatrix):
                        constraint = []
                        for j,(lbl,k) in enumerate(zip(coordVarLbl,row)):
                            if k == 0: continue
                            constraint.append([k,G2varObj[j]])
                        constraint += [shortmodelist[i],False,'f']
                        self.Constraints.append(constraint)
                    #----------------------------------------------------------------------
                    # save the ISODISTORT info for "mode analysis"
                    if 'ISODISTORT' not in self.Phase: self.Phase['ISODISTORT'] = {}
                    self.Phase['ISODISTORT'].update({
                        'IsoModeList' : modelist,
                        'G2ModeList' : shortmodelist,
                        'IsoVarList' : coordVarLbl,
                        'G2VarList' : G2varObj,
                        'ParentStructure' : ParentCoordinates,
                        'Var2ModeMatrix' : displacivemodeInvmatrix,
                        'Mode2VarMatrix' : displacivemodematrix,
                        })
                    # make explaination dictionary
                    for mode,shortmode in zip(modelist,shortmodelist):
                        explaination[shortmode] = "ISODISTORT full name "+str(mode)
                #----------------------------------------------------------------------
                # now read in the ISODISTORT occupancy modes
                #----------------------------------------------------------------------
                if blk.get('_iso_occupancymode_label'):
                    modelist = []
                    shortmodelist = []
                    for lbl in blk.get('_iso_occupancymode_label'):
                        modelist.append(lbl)
                        # assume lbl is of form SSSSS[x,y,z]AAAA(a,b,...)BBBBB
                        # where SSSSS is the parent spacegroup, [x,y,z] is a location 
                        regexp = re.match(r'.*?\[.*?\](.*?)\(.*?\)(.*)',lbl)
                        # this extracts the AAAAA and BBBBB parts of the string
                        if regexp:
                            lbl = regexp.expand(r'\1\2') # parse succeeded, make a short version
                        lbl = lbl.replace('order','o')
                        G2obj.MakeUniqueLabel(lbl,shortmodelist) # make unique and add to list
                    # read in the coordinate offset variables names and map them to G2 names/objects
                    occVarLbl = []
                    G2varLbl = []
                    G2varObj = []
                    error = False
                    for lbl in blk.get('_iso_deltaoccupancy_label'):
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
                        else:
                            anum = atomlbllist.index(albl)
                        var = self.varLookup.get(vlbl)
                        if not var:
                            self.warnings += ' ERROR: _iso_deltaoccupancy_label variable not found: '+lbl
                            error = True
                            continue
                        G2varLbl.append('::'+var+':'+str(anum)) # variable name, less phase ID
                        G2varObj.append(G2obj.G2VarObj(
                            (self.Phase['ranId'],None,var,ranIdlookup[albl])
                            ))
                    if error:
                        raise Exception,"Error decoding variable labels"

                    if len(G2varObj) != len(modelist):
                        print "non-square input"
                        raise Exception,"Rank of _iso_occupancymode != _iso_deltaoccupancy"

                    error = False
                    ParentCoordinates = {}
                    for lbl,exp in zip(
                        blk.get('_iso_occupancy_label'),
                        blk.get('_iso_occupancy_formula'),
                        ):
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
                            val = G2p3.FormulaEval(val)
                            if val is None:
                                self.warnings += ' ERROR: _iso_occupancy_formula coordinate not interpreted: '+lbl
                                error = True
                                continue
                            ParentCoordinates[albl] = val
                    if error:
                        raise Exception,"Error decoding occupancy labels"
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
                    for i,row in enumerate(occupancymodeInvmatrix):
                        constraint = []
                        for j,(lbl,k) in enumerate(zip(occVarLbl,row)):
                            if k == 0: continue
                            constraint.append([k,G2varObj[j]])
                        constraint += [shortmodelist[i],False,'f']
                        self.Constraints.append(constraint)
                    #----------------------------------------------------------------------
                    # save the ISODISTORT info for "mode analysis"
                    if 'ISODISTORT' not in self.Phase: self.Phase['ISODISTORT'] = {}
                    self.Phase['ISODISTORT'].update({
                        'OccModeList' : modelist,
                        'G2OccModeList' : shortmodelist,
                        'OccVarList' : occVarLbl,
                        'G2OccVarList' : G2varObj,
                        'BaseOcc' : ParentCoordinates,
                        'Var2OccMatrix' : occupancymodeInvmatrix,
                        'Occ2VarMatrix' : occupancymodematrix,
                        })
                    # make explaination dictionary
                    for mode,shortmode in zip(modelist,shortmodelist):
                        explaination[shortmode] = "ISODISTORT full name "+str(mode)
                #----------------------------------------------------------------------
                # done with read
                #----------------------------------------------------------------------
                if explaination: self.Constraints.append(explaination)

                # # debug: show the mode var to mode relations
                # for i,row in enumerate(displacivemodeInvmatrix):
                #     l = ''
                #     for j,(lbl,k) in enumerate(zip(coordVarLbl,row)):
                #         if k == 0: continue
                #         if l: l += ' + '
                #         #l += lbl+' * '+str(k)
                #         l += G2varLbl[j]+' * '+str(k)
                #     print str(i) + ': '+shortmodelist[i]+' = '+l
                # print 70*'='

                # # debug: Get the ISODISTORT offset values
                # coordVarDelta = {}
                # for lbl,val in zip(
                #     blk.get('_iso_deltacoordinate_label'),
                #     blk.get('_iso_deltacoordinate_value'),):
                #     coordVarDelta[lbl] = float(val)
                # modeVarDelta = {}
                # for lbl,val in zip(
                #     blk.get('_iso_displacivemode_label'),
                #     blk.get('_iso_displacivemode_value'),):
                #     modeVarDelta[lbl] = cif.get_number_with_esd(val)[0]

                # print 70*'='
                # # compute the mode values from the reported coordinate deltas
                # for i,row in enumerate(displacivemodeInvmatrix):
                #     l = ''
                #     sl = ''
                #     s = 0.
                #     for lbl,k in zip(coordVarLbl,row):
                #         if k == 0: continue
                #         if l: l += ' + '
                #         l += lbl+' * '+str(k)
                #         if sl: sl += ' + '
                #         sl += str(coordVarDelta[lbl])+' * '+str(k)
                #         s += coordVarDelta[lbl] * k
                #     print 'a'+str(i)+' = '+l
                #     print '\t= '+sl
                #     print  modelist[i],shortmodelist[i],modeVarDelta[modelist[i]],s
                #     print

                # print 70*'='
                # # compute the coordinate displacements from the reported mode values
                # for i,lbl,row in zip(range(len(coordVarLbl)),coordVarLbl,displacivemodematrix):
                #     l = ''
                #     sl = ''
                #     s = 0.0
                #     for j,k in enumerate(row):
                #         if k == 0: continue
                #         if l: l += ' + '
                #         l += 'a'+str(j+1)+' * '+str(k)
                #         if sl: sl += ' + '
                #         sl += str(shortmodelist[j]) +' = '+ str(modeVarDelta[modelist[j]]) + ' * '+str(k)
                #         s += modeVarDelta[modelist[j]] * k
                #     print lbl+' = '+l
                #     print '\t= '+sl
                #     print lbl,G2varLbl[i],coordVarDelta[lbl],s
                #     print

                # determine the coordinate delta values from deviations from the parent structure
                # for atmline in self.Phase['Atoms']:
                #     lbl = atmline[0]
                #     x,y,z = atmline[3:6]
                #     if lbl not in ParentCoordinates:
                #         print lbl,x,y,z
                #         continue
                #     px,py,pz = ParentCoordinates[lbl]
                #     print lbl,x,y,z,x-px,y-py,z-pz

                returnstat = True
        except Exception as detail:
            print 'CIF error:',detail # for testing
            print sys.exc_info()[0] # for testing
            import traceback
            print traceback.format_exc()
            returnstat = False
        finally:
            self.DoneBusy()
        return returnstat
