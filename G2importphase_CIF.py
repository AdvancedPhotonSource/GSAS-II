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
import GSASIIIO as G2IO
import GSASIIspc as G2spc
import GSASIIlattice as G2lat
import CifFile as cif # PyCifRW from James Hester
import urllib

class CIFPhaseReader(G2IO.ImportPhase):
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.CIF','.cif'),
            strictExtension=False,
            formatName = 'CIF',
            longFormatName = 'Crystallographic Information File import'
            )
    def ContentsValidator(self, filepointer):
        filepointer.seek(0) # rewind the file pointer
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
    def Reader(self,filename,filepointer, ParentFrame=None):
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
#### development code (to speed testing)
#            try:
#                    fp = open(filename+'cP',"rb")
#                    print("reading from "+filename+'cP')
#                    cf = cPickle.load(fp)
#                    fp.close()
#            except:
#                    cf = cif.ReadCif(filename)
#                    fp = open(filename+'cP',"wb")
#                    cPickle.dump(cf,fp)
#                    fp.close()
#### 
#### end development code
            self.ShowBusy() # this can take a while
            ciffile = 'file:'+urllib.pathname2url(filename)
            cf = cif.ReadCif(ciffile)
            #print cf
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
                for lbl in (
                    '_cell_length_a','_cell_length_b','_cell_length_c',
                    '_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma',
                    ):
                    cell.append(cif.get_number_with_esd(blk[lbl])[0])
                Volume = G2lat.calc_V(G2lat.cell2A(cell))
                self.Phase['General']['Cell'] = [False,]+cell+[Volume,]
                # read in atoms
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
                for aitem in atomloop:
                    atomlist = ['','','',0,0,0,1.0,'',0,'I',0.01,0,0,0,0,0,0]
                    atomlist.append(ran.randint(0,sys.maxint)) # add a unique label
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
