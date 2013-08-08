#!/usr/bin/env python
# -*- coding: utf-8 -*-
#G2cif
########### SVN repository information ###################
# $Date: 2013-07-22 20:57:37 -0500 (Mon, 22 Jul 2013) $
# $Author: toby $
# $Revision: 1006 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/exports/G2cif.py $
# $Id: G2cif.py 1006 2013-07-23 01:57:37Z toby $
########### SVN repository information ###################
'''Code to export a GSAS-II project as a CIF
The heavy lifting is done in method export
'''

# TODO: need a mechanism for editing of instrument names, bond pub flags, templates,...

import datetime as dt
import os.path
import numpy as np
import wx
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 1006 $")
import GSASIIIO as G2IO
#reload(G2IO)
import GSASIIgrid as G2gd
import GSASIIstrIO as G2stIO
#reload(G2stIO)
#import GSASIImapvars as G2mv
import GSASIImath as G2mth
#reload(G2mth)
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
#reload(G2spc)
import GSASIIphsGUI as G2pg
#reload(G2pg)
import GSASIIstrMain as G2stMn

DEBUG = False
#DEBUG = True    #True to skip printing of reflection/powder profile lists

def getCallerDocString(): # for development
    "Return the calling function's doc string"
    import inspect as ins
    for item in ins.stack()[1][0].f_code.co_consts:
        if type(item) is str:
            return item
    else:
        return '?'

class ExportCIF(G2IO.ExportBaseclass):
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'full CIF',
            extension='.cif',
            longFormatName = 'Export project as CIF'
            )
        self.author = ''

    def export(self,mode='full'):
        '''Export a CIF

        :param str mode: "full" (default) to create a complete CIF of project,
          "simple" for a simple CIF with only coordinates
        '''
    
        def openCIF(filnam):
            self.fp = open(filnam,'w')

        def closeCIF():
            self.fp.close()
            
        def WriteCIFitem(name,value=''):
            if value:
                if "\n" in value or len(value)> 70:
                    if name.strip(): self.fp.write(name+'\n')
                    self.fp.write('; '+value+'\n')
                    self.fp.write('; '+'\n')
                elif " " in value:
                    if len(name)+len(value) > 65:
                        self.fp.write(name + '\n   ' + '"' + str(value) + '"'+'\n')
                    else:
                        self.fp.write(name + '  ' + '"' + str(value) + '"'+'\n')
                else:
                    if len(name)+len(value) > 65:
                        self.fp.write(name+'\n   ' + value+'\n')
                    else:
                        self.fp.write(name+'  ' + value+'\n')
            else:
                self.fp.write(name+'\n')

        def WriteAudit():
            WriteCIFitem('_audit_creation_method',
                         'created in GSAS-II')
            WriteCIFitem('_audit_creation_date',self.CIFdate)
            if self.author:
                WriteCIFitem('_audit_author_name',self.author)
            WriteCIFitem('_audit_update_record',
                         self.CIFdate+'  Initial software-generated CIF')

        def WriteOverall():
            '''Write out overall refinement information.

            More could be done here, but this is a good start.
            '''
            WriteCIFitem('_pd_proc_info_datetime', self.CIFdate)
            WriteCIFitem('_pd_calc_method', 'Rietveld Refinement')
            #WriteCIFitem('_refine_ls_shift/su_max',DAT1)
            #WriteCIFitem('_refine_ls_shift/su_mean',DAT2)
            WriteCIFitem('_computing_structure_refinement','GSAS-II')
            try:
                vars = str(len(self.OverallParms['Covariance']['varyList']))
            except:
                vars = '?'
            WriteCIFitem('_refine_ls_number_parameters',vars)
            try:
                GOF = G2mth.ValEsd(self.OverallParms['Covariance']['Rvals']['GOF'],-0.009)
            except:
                GOF = '?'
            WriteCIFitem('_refine_ls_goodness_of_fit_all',GOF)

            # get restraint info
            # restraintDict = self.OverallParms.get('Restraints',{})
            # for i in  self.OverallParms['Constraints']:
            #     print i
            #     for j in self.OverallParms['Constraints'][i]:
            #         print j
            #WriteCIFitem('_refine_ls_number_restraints',TEXT)
            
            # other things to consider reporting
            # _refine_ls_number_reflns
            # _refine_ls_goodness_of_fit_obs
            # _refine_ls_R_factor_all
            # _refine_ls_R_factor_obs
            # _refine_ls_wR_factor_all
            # _refine_ls_wR_factor_obs
            # _refine_ls_restrained_S_all
            # _refine_ls_restrained_S_obs

            # include an overall profile r-factor, if there is more than one powder histogram
            if len(self.powderDict) > 1:
                WriteCIFitem('\n# OVERALL POWDER R-FACTOR')
                try:
                    R = str(self.OverallParms['Covariance']['Rvals']['Rwp'])
                except:
                    R = '?'
                WriteCIFitem('_pd_proc_ls_prof_wR_factor',R)
                #WriteCIFitem('_pd_proc_ls_prof_R_factor',TEXT(11:20)) # who cares!
            WriteCIFitem('_refine_ls_matrix_type','full')
            #WriteCIFitem('_refine_ls_matrix_type','userblocks')

        def WritePubTemplate():
            '''TODO: insert the publication template ``template_publ.cif`` or some modified
            version for this project. Store this in the GPX file?
            '''
            print getCallerDocString()

        def WritePhaseTemplate():
            '''TODO: insert the phase template ``template_phase.cif`` or some modified
            version for this project
            '''
            print getCallerDocString()

        def WritePowderTemplate():
            '''TODO: insert the phase template ``template_instrument.cif`` or some modified
            version for this project
            '''
            print getCallerDocString()

        def WriteSnglXtalTemplate():
            '''TODO: insert the single-crystal histogram template 
            for this project
            '''
            print getCallerDocString()

        def FormatSH(phasenam):
            'Format a full spherical harmonics texture description as a string'
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            pfx = str(phasedict['pId'])+'::'
            s = ""
            textureData = phasedict['General']['SH Texture']    
            if textureData.get('Order'):
                s += "Spherical Harmonics correction. Order = "+str(textureData['Order'])
                s += " Model: " + str(textureData['Model']) + "\n    Orientation angles: "
                for name in ['omega','chi','phi']:
                    aname = pfx+'SH '+name
                    s += name + " = "
                    sig = self.sigDict.get(aname,-0.09)
                    s += G2mth.ValEsd(self.parmDict[aname],sig)
                    s += "; "
                s += "\n"
                s1 = "    Coefficients:  "
                for name in textureData['SH Coeff'][1]:
                    aname = pfx+name
                    if len(s1) > 60:
                        s += s1 + "\n"
                        s1 = "    "
                    s1 += aname + ' = '
                    sig = self.sigDict.get(aname,-0.0009)
                    s1 += G2mth.ValEsd(self.parmDict[aname],sig)
                    s1 += "; "
                s += s1
            return s

        def FormatHAPpo(phasenam):
            '''return the March-Dollase/SH correction for every
            histogram in the current phase formatted into a
            character string
            '''
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            s = ''
            for histogram in sorted(phasedict['Histograms']):
                if histogram.startswith("HKLF"): continue # powder only
                Histogram = self.Histograms.get(histogram)
                if not Histogram: continue
                hapData = phasedict['Histograms'][histogram]
                if hapData['Pref.Ori.'][0] == 'MD':
                    aname = str(phasedict['pId'])+':'+str(Histogram['hId'])+':MD'
                    if self.parmDict.get(aname,1.0) != 1.0: continue
                    sig = self.sigDict.get(aname,-0.009)
                    if s != "": s += '\n'
                    s += 'March-Dollase correction'
                    if len(self.powderDict) > 1:
                        s += ', histogram '+str(Histogram['hId']+1)
                    s += ' coef. = ' + G2mth.ValEsd(self.parmDict[aname],sig)
                    s += ' axis = ' + str(hapData['Pref.Ori.'][3])
                else: # must be SH
                    if s != "": s += '\n'
                    s += 'Simple spherical harmonic correction'
                    if len(self.powderDict) > 1:
                        s += ', histogram '+str(Histogram['hId']+1)
                    s += ' Order = '+str(hapData['Pref.Ori.'][4])+'\n'
                    s1 = "    Coefficients:  "
                    for item in hapData['Pref.Ori.'][5]:
                        print item
                        aname = str(phasedict['pId'])+':'+str(Histogram['hId'])+':'+item
                        print aname
                        if len(s1) > 60:
                            s += s1 + "\n"
                            s1 = "    "
                        s1 += aname + ' = '
                        sig = self.sigDict.get(aname,-0.0009)
                        s1 += G2mth.ValEsd(self.parmDict[aname],sig)
                        s1 += "; "
                    s += s1
            return s
        def FormatBackground(bkg,hId):
            '''Display the Background information as a descriptive text string.
            
            TODO: this needs to be expanded to show the diffuse peak and
            Debye term information as well. (Bob)

            :returns: the text description (str)
            '''
            hfx = ':'+str(hId)+':'
            fxn, bkgdict = bkg
            terms = fxn[2]
            txt = 'Background function: "'+fxn[0]+'" function with '+str(terms)+' terms:\n'
            l = "    "
            for i,v in enumerate(fxn[3:]):
                name = '%sBack:%d'%(hfx,i)
                sig = self.sigDict.get(name,-0.009)
                if len(l) > 60:
                    txt += l + '\n'
                    l = '    '
                l += G2mth.ValEsd(v,sig)+', '
            txt += l
            if bkgdict['nDebye']:
                txt += '\n  Background Debye function parameters: A, R, U:'
                names = ['A:','R:','U:']
                for i in range(bkgdict['nDebye']):
                    txt += '\n    '
                    for j in range(3):
                        name = hfx+'Debye'+names[j]+str(i)
                        sig = self.sigDict.get(name,-0.009)
                        txt += G2mth.ValEsd(bkgdict['debyeTerms'][i][2*j],sig)+', '
            if bkgdict['nPeaks']:
                txt += '\n  Background peak parameters: pos, int, sig, gam:'
                names = ['pos:','int:','sig:','gam:']
                for i in range(bkgdict['nPeaks']):
                    txt += '\n    '
                    for j in range(4):
                        name = hfx+'BkPk'+names[j]+str(i)
                        sig = self.sigDict.get(name,-0.009)
                        txt += G2mth.ValEsd(bkgdict['peaksList'][i][2*j],sig)+', '
            return txt

        def FormatInstProfile(instparmdict,hId):
            '''Format the instrumental profile parameters with a
            string description. Will only be called on PWDR histograms
            '''
            s = ''
            inst = instparmdict[0]
            hfx = ':'+str(hId)+':'
            if 'C' in inst['Type'][0]:
                s = 'Finger-Cox-Jephcoat function parameters U, V, W, X, Y, SH/L:\n'
                s += '  peak variance(Gauss) = Utan(Th)^2+Vtan(Th)+W:\n'
                s += '  peak HW(Lorentz) = X/cos(Th)+Ytan(Th); SH/L = S/L+H/L\n'
                s += '  U, V, W in (centideg)^2, X & Y in centideg\n    '
                for item in ['U','V','W','X','Y','SH/L']:
                    name = hfx+item
                    sig = self.sigDict.get(name,-0.009)
                    s += G2mth.ValEsd(inst[item][1],sig)+', '                    
            elif 'T' in inst['Type'][0]:    #to be tested after TOF Rietveld done
                s = 'Von Dreele-Jorgenson-Windsor function parameters\n'+ \
                    '   alpha, beta-0, beta-1, beta-q, sig-0, sig-1, sig-q, X, Y:\n    '
                for item in ['alpha','bet-0','bet-1','bet-q','sig-0','sig-1','sig-q','X','Y']:
                    name = hfx+item
                    sig = self.sigDict.get(name,-0.009)
                    s += G2mth.ValEsd(inst[item][1],sig)+', '
            return s

        def FormatPhaseProfile(phasenam):
            '''Format the phase-related profile parameters (size/strain)
            with a string description.
            return an empty string or None if there are no
            powder histograms for this phase.
            '''
            s = ''
            phasedict = self.Phases[phasenam] # pointer to current phase info
            SGData = phasedict['General'] ['SGData']          
            for histogram in sorted(phasedict['Histograms']):
                if histogram.startswith("HKLF"): continue # powder only
                Histogram = self.Histograms.get(histogram)
                if not Histogram: continue
                hapData = phasedict['Histograms'][histogram]
                pId = phasedict['pId']
                hId = Histogram['hId']
                phfx = '%d:%d:'%(pId,hId)
                size = hapData['Size']
                mustrain = hapData['Mustrain']
                hstrain = hapData['HStrain']
                s = '  Crystallite size model "%s" for %s (microns)\n  '%(size[0],phasenam)
                names = ['Size;i','Size;mx']
                if 'uniax' in size[0]:
                    names = ['Size;i','Size;a','Size;mx']
                    s += 'anisotropic axis is %s\n  '%(str(size[3]))
                    s += 'parameters: equatorial size, axial size, G/L mix\n    '
                    for i,item in enumerate(names):
                        name = phfx+item
                        sig = self.sigDict.get(name,-0.009)
                        s += G2mth.ValEsd(size[1][i],sig)+', '
                elif 'ellip' in size[0]:
                    s += 'parameters: S11, S22, S33, S12, S13, S23, G/L mix\n    '
                    for i in range(6):
                        name = phfx+'Size:'+str(i)
                        sig = self.sigDict.get(name,-0.009)
                        s += G2mth.ValEsd(size[4][i],sig)+', '
                    sig = self.sigDict.get(phfx+'Size;mx',-0.009)
                    s += G2mth.ValEsd(size[1][2],sig)+', '                                           
                else:       #isotropic
                    s += 'parameters: Size, G/L mix\n    '
                    i = 0
                    for item in names:
                        name = phfx+item
                        sig = self.sigDict.get(name,-0.009)
                        s += G2mth.ValEsd(size[1][i],sig)+', '
                        i = 2    #skip the aniso value                
                s += '\n  Mustrain model "%s" for %s (10^6)\n  '%(mustrain[0],phasenam)
                names = ['Mustrain;i','Mustrain;mx']
                if 'uniax' in mustrain[0]:
                    names = ['Mustrain;i','Mustrain;a','Mustrain;mx']
                    s += 'anisotropic axis is %s\n  '%(str(size[3]))
                    s += 'parameters: equatorial mustrain, axial mustrain, G/L mix\n    '
                    for i,item in enumerate(names):
                        name = phfx+item
                        sig = self.sigDict.get(name,-0.009)
                        s += G2mth.ValEsd(mustrain[1][i],sig)+', '
                elif 'general' in mustrain[0]:
                    names = 'parameters: '
                    for i,name in enumerate(G2spc.MustrainNames(SGData)):
                        names += name+', '
                        if i == 9:
                            names += '\n  '
                    names += 'G/L mix\n    '
                    s += names
                    txt = ''
                    for i in range(len(mustrain[4])):
                        name = phfx+'Mustrain:'+str(i)
                        sig = self.sigDict.get(name,-0.009)
                        if len(txt) > 60:
                            s += txt+'\n    '
                            txt = ''
                        txt += G2mth.ValEsd(mustrain[4][i],sig)+', '
                    s += txt                                           
                    sig = self.sigDict.get(phfx+'Mustrain;mx',-0.009)
                    s += G2mth.ValEsd(mustrain[1][2],sig)+', '
                    
                else:       #isotropic
                    s += '  parameters: Mustrain, G/L mix\n    '
                    i = 0
                    for item in names:
                        name = phfx+item
                        sig = self.sigDict.get(name,-0.009)
                        s += G2mth.ValEsd(mustrain[1][i],sig)+', '
                        i = 2    #skip the aniso value                
                s += '\n  Macrostrain for %s\n'%(phasenam)
                txt = '  parameters: '
                names = G2spc.HStrainNames(SGData)
                for name in names:
                    txt += name+', '
                s += txt+'\n    '
                for i in range(len(names)):
                    name = phfx+name[i]
                    sig = self.sigDict.get(name,-0.009)
                    s += G2mth.ValEsd(hstrain[0][i],sig)+', '
            return s
        
        def FmtAtomType(sym):
            'Reformat a GSAS-II atom type symbol to match CIF rules'
            sym = sym.replace('_','') # underscores are not allowed: no isotope designation?
            # in CIF, oxidation state sign symbols come after, not before
            if '+' in sym:
                sym = sym.replace('+','') + '+'
            elif '-' in sym:
                sym = sym.replace('-','') + '-'
            return sym
            
        def PutInCol(val,wid):
            '''Pad a value to >=wid+1 columns by adding spaces at the end. Always
            adds at least one space
            '''
            val = str(val).replace(' ','')
            if not val: val = '?'
            fmt = '{:' + str(wid) + '} '
            return fmt.format(val)

        def MakeUniqueLabel(lbl,labellist):
            'Make sure that every atom label is unique'
            lbl = lbl.strip()
            if not lbl: # deal with a blank label
                lbl = 'A_1'
            if lbl not in labellist:
                labellist.append(lbl)
                return lbl
            i = 1
            prefix = lbl
            if '_' in lbl:
                prefix = lbl[:lbl.rfind('_')]
                suffix = lbl[lbl.rfind('_')+1:]
                try:
                    i = int(suffix)+1
                except:
                    pass
            while prefix+'_'+str(i) in labellist:
                i += 1
            else:
                lbl = prefix+'_'+str(i)
                labellist.append(lbl)

        def WriteAtomsNuclear(phasenam):
            'Write atom positions to CIF'
            phasedict = self.Phases[phasenam] # pointer to current phase info
            General = phasedict['General']
            cx,ct,cs,cia = General['AtomPtrs']
            Atoms = phasedict['Atoms']
            cfrac = cx+3
            fpfx = str(phasedict['pId'])+'::Afrac:'        
            for i,at in enumerate(Atoms):
                fval = self.parmDict.get(fpfx+str(i),at[cfrac])
                if fval != 0.0:
                    break
            else:
                WriteCIFitem('\n# PHASE HAS NO ATOMS!')
                return
                
            WriteCIFitem('\n# ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS')
            WriteCIFitem('loop_ '+
                         '\n\t_atom_site_label'+
                         '\n\t_atom_site_type_symbol'+
                         '\n\t_atom_site_fract_x'+
                         '\n\t_atom_site_fract_y'+
                         '\n\t_atom_site_fract_z'+
                         '\n\t_atom_site_occupancy'+
                         '\n\t_atom_site_adp_type'+
                         '\n\t_atom_site_U_iso_or_equiv'+
                         '\n\t_atom_site_symmetry_multiplicity')

            varnames = {cx:'Ax',cx+1:'Ay',cx+2:'Az',cx+3:'Afrac',
                        cia+1:'AUiso',cia+2:'AU11',cia+3:'AU22',cia+4:'AU33',
                        cia+5:'AU12',cia+6:'AU13',cia+7:'AU23'}
            self.labellist = []
            
            pfx = str(phasedict['pId'])+'::'
            # loop over all atoms
            naniso = 0
            for i,at in enumerate(Atoms):
                s = PutInCol(MakeUniqueLabel(at[ct-1],self.labellist),6) # label
                fval = self.parmDict.get(fpfx+str(i),at[cfrac])
                if fval == 0.0: continue # ignore any atoms that have a occupancy set to 0 (exact)
                s += PutInCol(FmtAtomType(at[ct]),4) # type
                if at[cia] == 'I':
                    adp = 'Uiso '
                else:
                    adp = 'Uani '
                    naniso += 1
                    # compute Uequiv crudely
                    # correct: Defined as "1/3 trace of diagonalized U matrix".
                    # SEE cell2GS & Uij2Ueqv to GSASIIlattice. Former is needed to make the GS matrix used by the latter.
                    t = 0.0
                    for j in (2,3,4):
                        var = pfx+varnames[cia+j]+":"+str(i)
                        t += self.parmDict.get(var,at[cia+j])
                for j in (cx,cx+1,cx+2,cx+3,cia+1):
                    if j in (cx,cx+1,cx+2):
                        dig = 11
                        sigdig = -0.00009
                    else:
                        dig = 10
                        sigdig = -0.009
                    var = pfx+varnames[j]+":"+str(i)
                    dvar = pfx+"d"+varnames[j]+":"+str(i)
                    if dvar not in self.sigDict:
                        dvar = var
                    if j == cia+1 and adp == 'Uani ':
                        val = t/3.
                        sig = sigdig
                    else:
                        #print var,(var in self.parmDict),(var in self.sigDict)
                        val = self.parmDict.get(var,at[j])
                        sig = self.sigDict.get(dvar,sigdig)
                    s += PutInCol(G2mth.ValEsd(val,sig),dig)
                s += adp
                s += PutInCol(at[cs+1],3)
                WriteCIFitem(s)
            if naniso == 0: return
            # now loop over aniso atoms
            WriteCIFitem('\nloop_' + '\n\t_atom_site_aniso_label' + 
                         '\n\t_atom_site_aniso_U_11' + '\n\t_atom_site_aniso_U_12' +
                         '\n\t_atom_site_aniso_U_13' + '\n\t_atom_site_aniso_U_22' +
                         '\n\t_atom_site_aniso_U_23' + '\n\t_atom_site_aniso_U_33')
            for i,at in enumerate(Atoms):
                fval = self.parmDict.get(fpfx+str(i),at[cfrac])
                if fval == 0.0: continue # ignore any atoms that have a occupancy set to 0 (exact)
                if at[cia] == 'I': continue
                s = PutInCol(self.labellist[i],6) # label
                for j in (2,3,4,5,6,7):
                    sigdig = -0.0009
                    var = pfx+varnames[cia+j]+":"+str(i)
                    val = self.parmDict.get(var,at[cia+j])
                    sig = self.sigDict.get(var,sigdig)
                    s += PutInCol(G2mth.ValEsd(val,sig),11)
                WriteCIFitem(s)

        def HillSortElements(elmlist):
            '''Sort elements in "Hill" order: C, H, others, (where others
            are alphabetical).

            :params list elmlist: a list of element strings

            :returns: a sorted list of element strings
            '''
            newlist = []
            oldlist = elmlist[:]
            for elm in ('C','H'):
                if elm in elmlist:
                    newlist.append(elm)
                    oldlist.pop(oldlist.index(elm))
            return newlist+sorted(oldlist)

        def WriteComposition(phasenam):
            '''determine the composition for the unit cell, crudely determine Z and
            then compute the composition in formula units
            '''
            phasedict = self.Phases[phasenam] # pointer to current phase info
            General = phasedict['General']
            Z = General.get('cellZ',0.0)
            cx,ct,cs,cia = General['AtomPtrs']
            Atoms = phasedict['Atoms']
            fpfx = str(phasedict['pId'])+'::Afrac:'        
            cfrac = cx+3
            cmult = cs+1
            compDict = {} # combines H,D & T
            sitemultlist = []
            massDict = dict(zip(General['AtomTypes'],General['AtomMass']))
            cellmass = 0
            for i,at in enumerate(Atoms):
                atype = at[ct].strip()
                if atype.find('-') != -1: atype = atype.split('-')[0]
                if atype.find('+') != -1: atype = atype.split('+')[0]
                atype = atype[0].upper()+atype[1:2].lower() # force case conversion
                if atype == "D" or atype == "D": atype = "H"
                fvar = fpfx+str(i)
                fval = self.parmDict.get(fvar,at[cfrac])
                mult = at[cmult]
                if not massDict.get(at[ct]):
                    print 'No mass found for atom type '+at[ct]
                    print 'Will not compute cell contents for phase '+phasenam
                    return
                cellmass += massDict[at[ct]]*mult*fval
                compDict[atype] = compDict.get(atype,0.0) + mult*fval
                if fval == 1: sitemultlist.append(mult)
            if len(compDict.keys()) == 0: return # no elements!
            if Z < 1: # Z has not been computed or set by user
                Z = 1
                for i in range(2,min(sitemultlist)+1):
                    for m in sitemultlist:
                        if m % i != 0:
                            break
                        else:
                            Z = i
                General['cellZ'] = Z # save it

            # when scattering factors are included in the CIF, this needs to be
            # added to the loop here but only in the one-block case.
            # For multiblock CIFs, scattering factors go in the histogram
            # blocks  (for all atoms in all appropriate phases)

            #if oneblock: # add scattering factors for current phase here
            WriteCIFitem('\nloop_  _atom_type_symbol _atom_type_number_in_cell')
            formula = ''
            reload(G2mth)
            for elem in HillSortElements(compDict.keys()):
                WriteCIFitem('  ' + PutInCol(elem,4) +
                             G2mth.ValEsd(compDict[elem],-0.009,True))
                if formula: formula += " "
                formula += elem
                if compDict[elem] == Z: continue
                formula += G2mth.ValEsd(compDict[elem]/Z,-0.009,True)
            WriteCIFitem( '\n# Note that Z affects _cell_formula_sum and _weight')
            WriteCIFitem( '_cell_formula_units_Z',str(Z))
            WriteCIFitem( '_chemical_formula_sum',formula)
            WriteCIFitem( '_chemical_formula_weight',
                          G2mth.ValEsd(cellmass/Z,-0.09,True))

        def WriteDistances(phasenam,SymOpList,offsetList,symOpList,G2oprList):
            '''Report bond distances and angles for the CIF

            Note that _geom_*_symmetry_* fields are values of form
            n_klm where n is the symmetry operation in SymOpList (counted
            starting with 1) and (k-5, l-5, m-5) are translations to add
            to (x,y,z). See
            http://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Igeom_angle_site_symmetry_.html

            TODO: need a method to select publication flags for distances/angles
            '''
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            Atoms = phasedict['Atoms']
            generalData = phasedict['General']
            cx,ct,cs,cia = phasedict['General']['AtomPtrs']
            cn = ct-1
            fpfx = str(phasedict['pId'])+'::Afrac:'        
            cfrac = cx+3
            DisAglData = {}
            DisAglCtls = {}
            xyz = []
            for i,atom in enumerate(Atoms):
                xyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
            if 'DisAglCtls' in generalData:
                DisAglCtls = generalData['DisAglCtls']
            else:
                dlg = G2gd.DisAglDialog(self.G2frame,DisAglCtls,generalData)
                if dlg.ShowModal() == wx.ID_OK:
                    DisAglCtls = dlg.GetData()
                    generalData['DisAglCtls'] = DisAglCtls
                else:
                    dlg.Destroy()
                    return
                dlg.Destroy()
            DisAglData['OrigAtoms'] = xyz
            DisAglData['TargAtoms'] = xyz
            DisAglData['SGData'] = generalData['SGData']
            SymOpList,offsetList,symOpList,G2oprList = G2spc.AllOps(
                generalData['SGData'])

#            print len(SymOpList),len(offsetList),len(symOpList),len(G2oprList)
#            raise Exception

            
            DisAglData['Cell'] = generalData['Cell'][1:] #+ volume
            if 'pId' in phasedict:
                DisAglData['pId'] = phasedict['pId']
                DisAglData['covData'] = self.OverallParms['Covariance']
            try:
                G2stMn.DistAngle(DisAglCtls,DisAglData)
            except KeyError:        # inside DistAngle for missing atom types in DisAglCtls
                print '**** ERROR - try again but do "Reset" to fill in missing atom types ****'
                    

            raise Exception
        
            # loop over interatomic distances for this phase
            WriteCIFitem('\n# MOLECULAR GEOMETRY')
            WriteCIFitem('loop_' + 
                         '\n\t_geom_bond_atom_site_label_1' +
                         '\n\t_geom_bond_atom_site_label_2' + 
                         '\n\t_geom_bond_distance' + 
                         '\n\t_geom_bond_site_symmetry_1' + 
                         '\n\t_geom_bond_site_symmetry_2' + 
                         '\n\t_geom_bond_publ_flag')

            # Note that labels should be read from self.labellist to correspond
            # to the values reported in the atoms table and zero occupancy atoms
            # should not be included
            fpfx = str(phasedict['pId'])+'::Afrac:'        
            for i,at in enumerate(Atoms):
                if self.parmDict.get(fpfx+str(i),at[cfrac]) == 0.0: continue
                lbl = self.labellist[i]


            # loop over interatomic angles for this phase
            WriteCIFitem('loop_' + 
                         '\n\t_geom_angle_atom_site_label_1' + 
                         '\n\t_geom_angle_atom_site_label_2' + 
                         '\n\t_geom_angle_atom_site_label_3' + 
                         '\n\t_geom_angle' + 
                         '\n\t_geom_angle_site_symmetry_1' +
                         '\n\t_geom_angle_site_symmetry_2' + 
                         '\n\t_geom_angle_site_symmetry_3' + 
                         '\n\t_geom_angle_publ_flag')


        def WritePhaseInfo(phasenam):
            WriteCIFitem('\n# phase info for '+str(phasenam) + ' follows')
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            WriteCIFitem('_pd_phase_name', phasenam)
            pfx = str(phasedict['pId'])+'::'
            A,sigA = G2stIO.cellFill(pfx,phasedict['General']['SGData'],self.parmDict,self.sigDict)
            cellSig = G2stIO.getCellEsd(pfx,
                                       phasedict['General']['SGData'],A,
                                       self.OverallParms['Covariance'])  # returns 7 vals, includes sigVol
            cellList = G2lat.A2cell(A) + (G2lat.calc_V(A),)
            defsigL = 3*[-0.00001] + 3*[-0.001] + [-0.01] # significance to use when no sigma
            names = ['length_a','length_b','length_c',
                     'angle_alpha','angle_beta ','angle_gamma',
                     'volume']
            prevsig = 0
            for lbl,defsig,val,sig in zip(names,defsigL,cellList,cellSig):
                if sig:
                    txt = G2mth.ValEsd(val,sig)
                    prevsig = -sig # use this as the significance for next value
                else:
                    txt = G2mth.ValEsd(val,min(defsig,prevsig),True)
                WriteCIFitem('_cell_'+lbl,txt)
                    
            WriteCIFitem('_symmetry_cell_setting',
                         phasedict['General']['SGData']['SGSys'])

            spacegroup = phasedict['General']['SGData']['SpGrp'].strip()
            # regularize capitalization and remove trailing H/R
            spacegroup = spacegroup[0].upper() + spacegroup[1:].lower().rstrip('rh ')
            WriteCIFitem('_symmetry_space_group_name_H-M',spacegroup)

            # generate symmetry operations including centering and center of symmetry
            SymOpList,offsetList,symOpList,G2oprList = G2spc.AllOps(
                phasedict['General']['SGData'])
            WriteCIFitem('loop_ _space_group_symop_id _space_group_symop_operation_xyz')
            for i,op in enumerate(SymOpList,start=1):
                WriteCIFitem('   {:3d}  {:}'.format(i,op.lower()))

            # loop over histogram(s) used in this phase
            if not oneblock and not self.quickmode:
                # report pointers to the histograms used in this phase
                histlist = []
                for hist in self.Phases[phasenam]['Histograms']:
                    if self.Phases[phasenam]['Histograms'][hist]['Use']:
                        if phasebyhistDict.get(hist):
                            phasebyhistDict[hist].append(phasenam)
                        else:
                            phasebyhistDict[hist] = [phasenam,]
                        blockid = datablockidDict.get(hist)
                        if not blockid:
                            print "Internal error: no block for data. Phase "+str(
                                phasenam)+" histogram "+str(hist)
                            histlist = []
                            break
                        histlist.append(blockid)

                if len(histlist) == 0:
                    WriteCIFitem('# Note: phase has no associated data')
                else:
                    WriteCIFitem('loop_  _pd_block_diffractogram_id')

            # report atom params
            if phasedict['General']['Type'] == 'nuclear':        #this needs macromolecular variant, etc!
                WriteAtomsNuclear(phasenam)
            else:
                raise Exception,"no export for mm coordinates implemented"
            # report cell contents
            WriteComposition(phasenam)
            if not self.quickmode:      # report distances and angles
                WriteDistances(phasenam,SymOpList,offsetList,symOpList,G2oprList)
                
        def Yfmt(ndec,val):
            'Format intensity values'
            out = ("{:."+str(ndec)+"f}").format(val)
            out = out.rstrip('0')  # strip zeros to right of decimal
            return out.rstrip('.')  # and decimal place when not needed

        def WritePowderData(histlbl):
            histblk = self.Histograms[histlbl]
            inst = histblk['Instrument Parameters'][0]
            hId = histblk['hId']
            pfx = ':' + str(hId) + ':'
            
            if 'Lam1' in inst:
                ratio = self.parmDict.get('I(L2)/I(L1)',inst['I(L2)/I(L1)'][1])
                sratio = self.sigDict.get('I(L2)/I(L1)',-0.0009)
                lam1 = self.parmDict.get('Lam1',inst['Lam1'][1])
                slam1 = self.sigDict.get('Lam1',-0.00009)
                lam2 = self.parmDict.get('Lam2',inst['Lam2'][1])
                slam2 = self.sigDict.get('Lam2',-0.00009)
                # always assume Ka1 & Ka2 if two wavelengths are present
                WriteCIFitem('loop_' + 
                             '\n\t_diffrn_radiation_wavelength' +
                             '\n\t_diffrn_radiation_wavelength_wt' + 
                             '\n\t_diffrn_radiation_type' + 
                             '\n\t_diffrn_radiation_wavelength_id')
                WriteCIFitem('  ' + PutInCol(G2mth.ValEsd(lam1,slam1),15)+
                             PutInCol('1.0',15) + 
                             PutInCol('K\\a~1~',10) + 
                             PutInCol('1',5))
                WriteCIFitem('  ' + PutInCol(G2mth.ValEsd(lam2,slam2),15)+
                             PutInCol(G2mth.ValEsd(ratio,sratio),15)+
                             PutInCol('K\\a~2~',10) + 
                             PutInCol('2',5))                
            else:
                lam1 = self.parmDict.get('Lam',inst['Lam'][1])
                slam1 = self.sigDict.get('Lam',-0.00009)
                WriteCIFitem('_diffrn_radiation_wavelength',G2mth.ValEsd(lam1,slam1))


            if not oneblock:
                if not phasebyhistDict.get(histlbl):
                    WriteCIFitem('\n# No phases associated with this data set')
                else:
                    WriteCIFitem('\n# PHASE TABLE')
                    WriteCIFitem('loop_' +
                                 '\n\t_pd_phase_id' + 
                                 '\n\t_pd_phase_block_id' + 
                                 '\n\t_pd_phase_mass_%')
                    wtFrSum = 0.
                    for phasenam in phasebyhistDict.get(histlbl):
                        hapData = self.Phases[phasenam]['Histograms'][histlbl]
                        General = self.Phases[phasenam]['General']
                        wtFrSum += hapData['Scale'][0]*General['Mass']

                    for phasenam in phasebyhistDict.get(histlbl):
                        hapData = self.Phases[phasenam]['Histograms'][histlbl]
                        General = self.Phases[phasenam]['General']
                        wtFr = hapData['Scale'][0]*General['Mass']/wtFrSum
                        pfx = str(self.Phases[phasenam]['pId'])+':'+str(hId)+':'
                        if pfx+'Scale' in self.sigDict:
                            sig = self.sigDict[pfx+'Scale']*wtFr/hapData['Scale'][0]
                        else:
                            sig = -0.0001
                        WriteCIFitem(
                            '  '+
                            str(self.Phases[phasenam]['pId']) +
                            '  '+datablockidDict[phasenam]+
                            '  '+G2mth.ValEsd(wtFr,sig)
                            )

            # TODO: this will need help from Bob
            # WriteCIFitem('_pd_proc_ls_prof_R_factor','?')
            # WriteCIFitem('_pd_proc_ls_prof_wR_factor','?')
            # WriteCIFitem('_pd_proc_ls_prof_wR_expected','?')
            # WriteCIFitem('_refine_ls_R_Fsqd_factor','?')

            if histblk['Instrument Parameters'][0]['Type'][1][1] == 'X':
                WriteCIFitem('_diffrn_radiation_probe','x-ray')
                pola = histblk['Instrument Parameters'][0].get('Polariz.')
                if pola:
                    pfx = ':' + str(hId) + ':'
                    sig = self.sigDict.get(pfx+'Polariz.',-0.0009)
                    txt = G2mth.ValEsd(pola[1],sig)
                    WriteCIFitem('_diffrn_radiation_polarisn_ratio',txt)
            elif histblk['Instrument Parameters'][0]['Type'][1][1] == 'N':
                WriteCIFitem('_diffrn_radiation_probe','neutron')
            # TOF (note that this may not be defined)
            #if histblk['Instrument Parameters'][0]['Type'][1][2] == 'T':
            #    WriteCIFitem('_pd_meas_2theta_fixed',text)
            

            # TODO: this will need help from Bob
            #if not oneblock:
            #WriteCIFitem('\n# SCATTERING FACTOR INFO')
            #WriteCIFitem('loop_  _atom_type_symbol')
            #if histblk['Instrument Parameters'][0]['Type'][1][1] == 'X':
            #    WriteCIFitem('      _atom_type_scat_dispersion_real')
            #    WriteCIFitem('      _atom_type_scat_dispersion_imag')
            #    for lbl in ('a1','a2','a3', 'a4', 'b1', 'b2', 'b3', 'b4', 'c'):
            #        WriteCIFitem('      _atom_type_scat_Cromer_Mann_'+lbl)
            #elif histblk['Instrument Parameters'][0]['Type'][1][1] == 'N':
            #    WriteCIFitem('      _atom_type_scat_length_neutron')
            #WriteCIFitem('      _atom_type_scat_source')

            WriteCIFitem('_pd_proc_ls_background_function',FormatBackground(histblk['Background'],histblk['hId']))

            # TODO: this will need help from Bob
            #WriteCIFitem('_exptl_absorpt_process_details','?')
            #WriteCIFitem('_exptl_absorpt_correction_T_min','?')
            #WriteCIFitem('_exptl_absorpt_correction_T_max','?')
            #C extinction
            #WRITE(IUCIF,'(A)') '# Extinction correction'
            #CALL WRVAL(IUCIF,'_gsas_exptl_extinct_corr_T_min',TEXT(1:10))
            #CALL WRVAL(IUCIF,'_gsas_exptl_extinct_corr_T_max',TEXT(11:20))

            if not oneblock:                 # instrumental profile terms go here
                WriteCIFitem('_pd_proc_ls_profile_function', 
                    FormatInstProfile(histblk["Instrument Parameters"],histblk['hId']))

            #refprx = '_refln.' # mm
            refprx = '_refln_' # normal
            WriteCIFitem('\n# STRUCTURE FACTOR TABLE')            
            # compute maximum intensity reflection
            Imax = 0
            for phasenam in histblk['Reflection Lists']:
                scale = self.Phases[phasenam]['Histograms'][histlbl]['Scale'][0]
                Icorr = np.array([refl[13] for refl in histblk['Reflection Lists'][phasenam]])
                FO2 = np.array([refl[8] for refl in histblk['Reflection Lists'][phasenam]])
                I100 = scale*FO2*Icorr
                Imax = max(Imax,max(I100))

            WriteCIFitem('loop_')
            if len(histblk['Reflection Lists'].keys()) > 1:
                WriteCIFitem('\t_pd_refln_phase_id')
            WriteCIFitem('\t' + refprx + 'index_h' + 
                         '\n\t' + refprx + 'index_k' + 
                         '\n\t' + refprx + 'index_l' + 
                         '\n\t' + refprx + 'F_squared_meas' + 
                         '\n\t' + refprx + 'F_squared_calc' + 
                         '\n\t' + refprx + 'phase_calc' + 
                         '\n\t_pd_refln_d_spacing')
            if Imax > 0:
                WriteCIFitem('\t_gsas_i100_meas')

            refcount = 0
            hklmin = None
            hklmax = None
            dmax = None
            dmin = None
            for phasenam in histblk['Reflection Lists']:
                scale = self.Phases[phasenam]['Histograms'][histlbl]['Scale'][0]
                phaseid = self.Phases[phasenam]['pId']
                refcount += len(histblk['Reflection Lists'][phasenam])
                for ref in histblk['Reflection Lists'][phasenam]:
                    if DEBUG:
                        print 'DEBUG: skip reflection list'
                        break
                    if hklmin is None:
                        hklmin = ref[0:3]
                        hklmax = ref[0:3]
                        dmax = dmin = ref[4]
                    if len(histblk['Reflection Lists'].keys()) > 1:
                        s = PutInCol(phaseid,2)
                    else:
                        s = ""
                    for i,hkl in enumerate(ref[0:3]):
                        hklmax[i] = max(hkl,hklmax[i])
                        hklmin[i] = min(hkl,hklmin[i])
                        s += PutInCol(int(hkl),4)
                    for I in ref[8:10]:
                        s += PutInCol(G2mth.ValEsd(I,-0.0009),10)
                    s += PutInCol(G2mth.ValEsd(ref[10],-0.9),7)
                    dmax = max(dmax,ref[4])
                    dmin = min(dmin,ref[4])
                    s += PutInCol(G2mth.ValEsd(ref[4],-0.009),8)
                    if Imax > 0:
                        I100 = 100.*scale*ref[8]*ref[13]/Imax
                        s += PutInCol(G2mth.ValEsd(I100,-0.09),6)
                    WriteCIFitem("  "+s)

            WriteCIFitem('_reflns_number_total', str(refcount))
            if hklmin is not None and len(histblk['Reflection Lists']) == 1: # hkl range has no meaning with multiple phases
                WriteCIFitem('_reflns_limit_h_min', str(int(hklmin[0])))
                WriteCIFitem('_reflns_limit_h_max', str(int(hklmax[0])))
                WriteCIFitem('_reflns_limit_k_min', str(int(hklmin[1])))
                WriteCIFitem('_reflns_limit_k_max', str(int(hklmax[1])))
                WriteCIFitem('_reflns_limit_l_min', str(int(hklmin[2])))
                WriteCIFitem('_reflns_limit_l_max', str(int(hklmax[2])))
            if hklmin is not None:
                WriteCIFitem('_reflns_d_resolution_high', G2mth.ValEsd(dmin,-0.009))
                WriteCIFitem('_reflns_d_resolution_low', G2mth.ValEsd(dmax,-0.0009))

            WriteCIFitem('\n# POWDER DATA TABLE')
            # is data fixed step? If the step varies by <0.01% treat as fixed step
            steps = histblk['Data'][0][1:] - histblk['Data'][0][:-1]
            if abs(max(steps)-min(steps)) > abs(max(steps))/10000.:
                fixedstep = False
            else:
                fixedstep = True

            if fixedstep: # and not TOF
                WriteCIFitem('_pd_meas_2theta_range_min', G2mth.ValEsd(histblk['Data'][0][0],-0.00009))
                WriteCIFitem('_pd_meas_2theta_range_max', G2mth.ValEsd(histblk['Data'][0][-1],-0.00009))
                WriteCIFitem('_pd_meas_2theta_range_inc', G2mth.ValEsd(steps.sum()/len(steps),-0.00009))
                # zero correct, if defined
                zero = None
                zerolst = histblk['Instrument Parameters'][0].get('Zero')
                if zerolst: zero = zerolst[1]
                zero = self.parmDict.get('Zero',zero)
                if zero:
                    WriteCIFitem('_pd_proc_2theta_range_min', G2mth.ValEsd(histblk['Data'][0][0]-zero,-0.00009))
                    WriteCIFitem('_pd_proc_2theta_range_max', G2mth.ValEsd(histblk['Data'][0][-1]-zero,-0.00009))
                    WriteCIFitem('_pd_proc_2theta_range_inc', G2mth.ValEsd(steps.sum()/len(steps),-0.00009))
                
            if zero:
                WriteCIFitem('_pd_proc_number_of_points', str(len(histblk['Data'][0])))
            else:
                WriteCIFitem('_pd_meas_number_of_points', str(len(histblk['Data'][0])))
            WriteCIFitem('\nloop_')
            #            WriteCIFitem('\t_pd_proc_d_spacing') # need easy way to get this
            if not fixedstep:
                if zero:
                    WriteCIFitem('\t_pd_proc_2theta_corrected')
                else:
                    WriteCIFitem('\t_pd_meas_2theta_scan')
            # at least for now, always report weights.
            #if countsdata:
            #    WriteCIFitem('\t_pd_meas_counts_total')
            #else:
            WriteCIFitem('\t_pd_meas_intensity_total')
            WriteCIFitem('\t_pd_calc_intensity_total')
            WriteCIFitem('\t_pd_proc_intensity_bkg_calc')
            WriteCIFitem('\t_pd_proc_ls_weight')
            maxY = max(histblk['Data'][1].max(),histblk['Data'][3].max())
            if maxY < 0: maxY *= -10 # this should never happen, but...
            ndec = max(0,10-int(np.log10(maxY))-1) # 10 sig figs should be enough
            maxSU = histblk['Data'][2].max()
            if maxSU < 0: maxSU *= -1 # this should never happen, but...
            ndecSU = max(0,8-int(np.log10(maxSU))-1) # 8 sig figs should be enough
            lowlim,highlim = histblk['Limits'][1]

            for x,yobs,yw,ycalc,ybkg in zip(histblk['Data'][0],
                                            histblk['Data'][1],
                                            histblk['Data'][2],
                                            histblk['Data'][3],
                                            histblk['Data'][4]):
                if lowlim <= x <= highlim:
                    pass
                else:
                    yw = 0.0 # show the point is not in use

                if fixedstep:
                    s = ""
                else:
                    s = PutInCol(G2mth.ValEsd(x-zero,-0.00009),10)
                s += PutInCol(Yfmt(ndec,yobs),12)
                s += PutInCol(Yfmt(ndec,ycalc),12)
                s += PutInCol(Yfmt(ndec,ybkg),11)
                s += PutInCol(Yfmt(ndecSU,yw),9)
                WriteCIFitem("  "+s)

        def WriteSingleXtalData(histlbl):
            histblk = self.Histograms[histlbl]
            #refprx = '_refln.' # mm
            refprx = '_refln_' # normal

            WriteCIFitem('\n# STRUCTURE FACTOR TABLE')            
            WriteCIFitem('loop_' + 
                         '\n\t' + refprx + 'index_h' + 
                         '\n\t' + refprx + 'index_k' + 
                         '\n\t' + refprx + 'index_l' +
                         '\n\t' + refprx + 'F_squared_meas' + 
                         '\n\t' + refprx + 'F_squared_sigma' + 
                         '\n\t' + refprx + 'F_squared_calc' + 
                         '\n\t' + refprx + 'phase_calc' +
                         '\n\t_pd_refln_d_spacing'
                         )

            hklmin = None
            hklmax = None
            dmax = None
            dmin = None
            refcount = len(histblk['Data'])
            for ref in histblk['Data']:
                s = "  "
                if hklmin is None:
                    hklmin = ref[0:3]
                    hklmax = ref[0:3]
                    dmax = dmin = ref[4]
                for i,hkl in enumerate(ref[0:3]):
                    hklmax[i] = max(hkl,hklmax[i])
                    hklmin[i] = min(hkl,hklmin[i])
                    s += PutInCol(int(hkl),4)
                sig = ref[6] * ref[8] / ref[5]
                s += PutInCol(G2mth.ValEsd(ref[8],-abs(sig/10)),12)
                s += PutInCol(G2mth.ValEsd(sig,-abs(sig)/10.),10)
                s += PutInCol(G2mth.ValEsd(ref[9],-abs(sig/10)),12)
                s += PutInCol(G2mth.ValEsd(ref[10],-0.9),7)
                dmax = max(dmax,ref[4])
                dmin = min(dmin,ref[4])
                s += PutInCol(G2mth.ValEsd(ref[4],-0.009),8)
                WriteCIFitem(s)
            WriteCIFitem('_reflns_number_total', str(refcount))
            if hklmin is not None:
                WriteCIFitem('_reflns_limit_h_min', str(int(hklmin[0])))
                WriteCIFitem('_reflns_limit_h_max', str(int(hklmax[0])))
                WriteCIFitem('_reflns_limit_k_min', str(int(hklmin[1])))
                WriteCIFitem('_reflns_limit_k_max', str(int(hklmax[1])))
                WriteCIFitem('_reflns_limit_l_min', str(int(hklmin[2])))
                WriteCIFitem('_reflns_limit_l_max', str(int(hklmax[2])))
                WriteCIFitem('_reflns_d_resolution_high', G2mth.ValEsd(dmin,-0.009))
                WriteCIFitem('_reflns_d_resolution_low', G2mth.ValEsd(dmax,-0.0009))

        #============================================================
        # the export process starts here
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()

        # Someday: get restraint & constraint info
        #restraintDict = self.OverallParms.get('Restraints',{})
        #for i in  self.OverallParms['Constraints']:
        #    print i
        #    for j in self.OverallParms['Constraints'][i]:
        #        print j

        self.CIFdate = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%dT%H:%M")
        # index powder and single crystal histograms
        self.powderDict = {}
        self.xtalDict = {}
        for hist in self.Histograms:
            i = self.Histograms[hist]['hId']
            if hist.startswith("PWDR"): 
                self.powderDict[i] = hist
            elif hist.startswith("HKLF"): 
                self.xtalDict[i] = hist
        # is there anything to export?
        if len(self.Phases) + len(self.powderDict) + len(self.xtalDict) == 0: 
            self.G2frame.ErrorDialog(
                'Empty project',
                'No data or phases to include in CIF')
            return
        # is there a file name defined?
        self.CIFname = os.path.splitext(
            os.path.split(self.G2frame.GSASprojectfile)[1]
            )[0]
        self.CIFname = self.CIFname.replace(' ','')
        if not self.CIFname:
            self.G2frame.ErrorDialog(
                'No GPX name',
                'Please save the project to provide a name')
            return
        # test for quick CIF mode or no data
        self.quickmode = False
        phasenam = phasenum = None # include all phases
        if mode != "full" or len(self.powderDict) + len(self.xtalDict) == 0:
            self.quickmode = True
            oneblock = True
            if len(self.Phases) == 0:
                self.G2frame.ErrorDialog(
                    'No phase present',
                    'Cannot create a coordinates CIF with no phases')
                return
            elif len(self.Phases) > 1: # quick mode: choose one phase
                choices = sorted(self.Phases.keys())
                phasenum = G2gd.ItemSelector(choices,self.G2frame)
                if phasenum is None: return
                phasenam = choices[phasenum]
        # will this require a multiblock CIF?
        elif len(self.Phases) > 1:
            oneblock = False
        elif len(self.powderDict) + len(self.xtalDict) > 1:
            oneblock = False
        else: # one phase, one dataset, Full CIF
            oneblock = True

        # make sure needed infomation is present
        # get CIF author name -- required for full CIFs
        try:
            self.author = self.OverallParms['Controls'].get("Author",'').strip()
        except KeyError:
            pass
        while not (self.author or self.quickmode):
            dlg = G2gd.SingleStringDialog(self.G2frame,'Get CIF Author','Provide CIF Author name (Last, First)')
            if not dlg.Show(): return # cancel was pressed
            self.author = dlg.GetValue()
            dlg.Destroy()
        try:
            self.OverallParms['Controls']["Author"] = self.author # save for future
        except KeyError:
            pass
        self.shortauthorname = self.author.replace(',','').replace(' ','')[:20]

        # check the instrument name for every histogram
        if not self.quickmode:
            dictlist = []
            keylist = []
            lbllist = []
            invalid = 0
            key3 = 'InstrName'
            for hist in self.Histograms:
                if hist.startswith("PWDR"): 
                    key2 = "Sample Parameters"
                    d = self.Histograms[hist][key2]
                elif hist.startswith("HKLF"): 
                    key2 = "Instrument Parameters"
                    d = self.Histograms[hist][key2][0]
                    
                lbllist.append(hist)
                dictlist.append(d)
                keylist.append(key3)
                instrname = d.get(key3)
                if instrname is None:
                    d[key3] = ''
                    invalid += 1
                elif instrname.strip() == '':
                    invalid += 1
            if invalid:
                msg = ""
                if invalid > 3: msg = (
                    "\n\nNote: it may be faster to set the name for\n"
                    "one histogram for each instrument and use the\n"
                    "File/Copy option to duplicate the name"
                    )
                if not G2gd.CallScrolledMultiEditor(
                    self.G2frame,dictlist,keylist,
                    prelbl=range(1,len(dictlist)+1),
                    postlbl=lbllist,
                    title='Instrument names',
                    header="Edit instrument names. Note that a non-blank\nname is required for all histograms"+msg,
                    ): return

        if oneblock and not self.quickmode:
            # select a dataset to use (there should only be one set in one block,
            # but take whatever comes 1st
            for hist in self.Histograms:
                histblk = self.Histograms[hist]
                if hist.startswith("PWDR"): 
                    instnam = histblk["Sample Parameters"]['InstrName']
                    break # ignore all but 1st data histogram
                elif hist.startswith("HKLF"): 
                    instnam = histblk["Instrument Parameters"][0]['InstrName']
                    break # ignore all but 1st data histogram
        if self.quickmode:
            fil = self.askSaveFile()
        else:
            fil = self.defSaveFile()
        if not fil: return
        openCIF(fil)
        #======================================================================
        # Start writing the CIF - single block
        #======================================================================
        if oneblock:
            WriteCIFitem('data_'+self.CIFname)
            if phasenam is None: # if not already selected, select the first phase (should be one) 
                phasenam = self.Phases.keys()[0]
            #print 'phasenam',phasenam
            phaseblk = self.Phases[phasenam] # pointer to current phase info
            if not self.quickmode:
                instnam = instnam.replace(' ','')
                WriteCIFitem('_pd_block_id',
                             str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                             str(self.shortauthorname) + "|" + instnam)
                WriteAudit()
                WritePubTemplate()
                WriteOverall()
                WritePhaseTemplate()
            # report the phase info
            WritePhaseInfo(phasenam)
            if hist.startswith("PWDR") and not self.quickmode:
                # preferred orientation
                SH = FormatSH(phasenam)
                MD = FormatHAPpo(phasenam)
                if SH and MD:
                    WriteCIFitem('_pd_proc_ls_pref_orient_corr', SH + '\n' + MD)
                elif SH or MD:
                    WriteCIFitem('_pd_proc_ls_pref_orient_corr', SH + MD)
                else:
                    WriteCIFitem('_pd_proc_ls_pref_orient_corr', 'none')
                    # report profile, since one-block: include both histogram and phase info
                WriteCIFitem('_pd_proc_ls_profile_function',
                    FormatInstProfile(histblk["Instrument Parameters"],histblk['hId'])
                    +'\n'+FormatPhaseProfile(phasenam))
                WritePowderTemplate()
                WritePowderData(hist)
            elif hist.startswith("HKLF") and not self.quickmode:
                WriteSnglXtalTemplate()
                WriteSingleXtalData(hist)
        else:
        #======================================================================
        # Start writing the CIF - multiblock
        #======================================================================
            # publication info
            WriteCIFitem('\ndata_'+self.CIFname+'_publ')
            WriteAudit()
            WriteCIFitem('_pd_block_id',
                         str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                         str(self.shortauthorname) + "|Overall")
            WritePubTemplate()
            # overall info
            WriteCIFitem('data_'+str(self.CIFname)+'_overall')
            WriteOverall()
            #============================================================
            WriteCIFitem('# POINTERS TO PHASE AND HISTOGRAM BLOCKS')
            datablockidDict = {} # save block names here -- N.B. check for conflicts between phase & hist names (unlikely!)
            # loop over phase blocks
            if len(self.Phases) > 1:
                loopprefix = ''
                WriteCIFitem('loop_   _pd_phase_block_id')
            else:
                loopprefix = '_pd_phase_block_id'
            
            for phasenam in sorted(self.Phases.keys()):
                i = self.Phases[phasenam]['pId']
                datablockidDict[phasenam] = (str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                             'phase_'+ str(i) + '|' + str(self.shortauthorname))
                WriteCIFitem(loopprefix,datablockidDict[phasenam])
            # loop over data blocks
            if len(self.powderDict) + len(self.xtalDict) > 1:
                loopprefix = ''
                WriteCIFitem('loop_   _pd_block_diffractogram_id')
            else:
                loopprefix = '_pd_block_diffractogram_id'
            for i in sorted(self.powderDict.keys()):
                hist = self.powderDict[i]
                histblk = self.Histograms[hist]
                instnam = histblk["Sample Parameters"]['InstrName']
                instnam = instnam.replace(' ','')
                i = histblk['hId']
                datablockidDict[hist] = (str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                                         str(self.shortauthorname) + "|" +
                                         instnam + "_hist_"+str(i))
                WriteCIFitem(loopprefix,datablockidDict[hist])
            for i in sorted(self.xtalDict.keys()):
                hist = self.xtalDict[i]
                histblk = self.Histograms[hist]
                instnam = histblk["Instrument Parameters"][0]['InstrName']
                instnam = instnam.replace(' ','')
                i = histblk['hId']
                datablockidDict[hist] = (str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                                         str(self.shortauthorname) + "|" +
                                         instnam + "_hist_"+str(i))
                WriteCIFitem(loopprefix,datablockidDict[hist])
            #============================================================
            # loop over phases, exporting them
            phasebyhistDict = {} # create a cross-reference to phases by histogram
            for j,phasenam in enumerate(sorted(self.Phases.keys())):
                i = self.Phases[phasenam]['pId']
                WriteCIFitem('\ndata_'+self.CIFname+"_phase_"+str(i))
                print "debug, processing ",phasenam
                WriteCIFitem('# Information for phase '+str(i))
                WriteCIFitem('_pd_block_id',datablockidDict[phasenam])
                # report the phase
                WritePhaseTemplate()
                WritePhaseInfo(phasenam)
                # preferred orientation
                SH = FormatSH(phasenam)
                MD = FormatHAPpo(phasenam)
                if SH and MD:
                    WriteCIFitem('_pd_proc_ls_pref_orient_corr', SH + '\n' + MD)
                elif SH or MD:
                    WriteCIFitem('_pd_proc_ls_pref_orient_corr', SH + MD)
                else:
                    WriteCIFitem('_pd_proc_ls_pref_orient_corr', 'none')
                # report sample profile terms
                PP = FormatPhaseProfile(phasenam)
                if PP:
                    WriteCIFitem('_pd_proc_ls_profile_function',PP)
                    
            #============================================================
            # loop over histograms, exporting them
            for i in sorted(self.powderDict.keys()):
                hist = self.powderDict[i]
                histblk = self.Histograms[hist]
                if hist.startswith("PWDR"): 
                    WriteCIFitem('\ndata_'+self.CIFname+"_pwd_"+str(i))
                    #instnam = histblk["Sample Parameters"]['InstrName']
                    # report instrumental profile terms
                    WriteCIFitem('_pd_proc_ls_profile_function',
                        FormatInstProfile(histblk["Instrument Parameters"],histblk['hId']))
                    WriteCIFitem('# Information for histogram '+str(i)+': '+hist)
                    WriteCIFitem('_pd_block_id',datablockidDict[hist])
                    WritePowderTemplate()
                    WritePowderData(hist)
            for i in sorted(self.xtalDict.keys()):
                hist = self.xtalDict[i]
                histblk = self.Histograms[hist]
                if hist.startswith("HKLF"): 
                    WriteCIFitem('\ndata_'+self.CIFname+"_sx_"+str(i))
                    #instnam = histblk["Instrument Parameters"][0]['InstrName']
                    WriteCIFitem('# Information for histogram '+str(i)+': '+hist)
                    WriteCIFitem('_pd_block_id',datablockidDict[hist])
                    WriteSnglXtalTemplate()
                    WriteSingleXtalData(hist)

        WriteCIFitem('#--' + 15*'eof--' + '#')
        closeCIF()
