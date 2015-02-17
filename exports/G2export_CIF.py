#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module G2export_CIF: CIF Exports*
------------------------------------------------------

This implements a complex exporter :class:`ExportCIF` that can implement an
entire project in a complete CIF intended for submission as a
publication. In addition, there are two subclasses of :class:`ExportCIF`:
:class:`ExportPhaseCIF` and :class:`ExportDataCIF` that
export a single phase or data set. Note that ``self.mode`` determines
what is written:

 * `self.mode="simple"` creates a simple CIF with only coordinates
   or data, while 

 * `self.mode="full"` creates a complete CIF of project.
 
'''

import datetime as dt
import os.path
import sys
import numpy as np
import cPickle
import copy
import re
import wx
import wx.lib.scrolledpanel as wxscroll
import wx.lib.resizewidget as rw
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIIO as G2IO
import GSASIIgrid as G2gd
import GSASIIctrls as G2ctrls
import GSASIIstrIO as G2stIO
import GSASIImath as G2mth
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIphsGUI as G2pg
import GSASIIstrMain as G2stMn

DEBUG = False    #True to skip printing of reflection/powder profile lists

CIFdic = None

class ExportCIF(G2IO.ExportBaseclass):
    '''Used to create a CIF of an entire project

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'Full CIF',
            extension='.cif',
            longFormatName = 'Export project as CIF'
            )
        self.exporttype = ['project']
        self.author = ''
        self.mode = 'full'

    def Exporter(self,event=None):
        '''Export a CIF. Export can be full or simple (as set by self.mode).
        "simple" skips data, distances & angles, etc. and can only include
        a single phase while "full" is intended for for publication submission.
        '''

#***** define functions for export method =======================================
        def WriteCIFitem(name,value=''):
            '''Write CIF data items to the file. Formats values as needed.
            Also used without a value for loops, comments, loop headers, etc.
            '''            
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
            'Write the CIF audit values. Perhaps should be in a single element loop.'
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
            if self.ifPWDR:
                WriteCIFitem('_pd_proc_info_datetime', self.CIFdate)
                WriteCIFitem('_pd_calc_method', 'Rietveld Refinement')
            #WriteCIFitem('_refine_ls_shift/su_max',DAT1)
            #WriteCIFitem('_refine_ls_shift/su_mean',DAT2)
            WriteCIFitem('_computing_structure_refinement','GSAS-II (Toby & Von Dreele, J. Appl. Cryst. 46, 544-549, 2013)')
            if self.ifHKLF:
                controls = self.OverallParms['Controls']
                if controls['F**2']:
                    thresh = 'F**2>%.1fu(F**2)'%(controls['minF/sig'])
                else:
                    thresh = 'F>%.1fu(F)'%(controls['minF/sig'])
                WriteCIFitem('_reflns_threshold_expression', thresh)                
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
            # _refine_ls_wR_factor_obs
            # _refine_ls_restrained_S_all
            # _refine_ls_restrained_S_obs

            # include an overall profile r-factor, if there is more than one powder histogram
            R = '%.5f'%(self.OverallParms['Covariance']['Rvals']['Rwp']/100.)
            WriteCIFitem('\n# OVERALL WEIGHTED R-FACTOR')
            WriteCIFitem('_refine_ls_wR_factor_obs',R)
                # _refine_ls_R_factor_all
                # _refine_ls_R_factor_obs                
            WriteCIFitem('_refine_ls_matrix_type','full')
            #WriteCIFitem('_refine_ls_matrix_type','userblocks')

        def writeCIFtemplate(G2dict,tmplate,defaultname=''):
            '''Write out the selected or edited CIF template
            An unedited CIF template file is copied, comments intact; an edited
            CIF template is written out from PyCifRW which of course strips comments. 
            In all cases the initial data_ header is stripped (there should only be one!)
            '''
            CIFobj = G2dict.get("CIF_template")
            if defaultname:
                defaultname = defaultname.encode('ascii','replace').strip().replace(' ','_')
                defaultname = re.sub(r'[^a-zA-Z0-9_-]','',defaultname)
                defaultname = tmplate + "_" + defaultname + ".cif"
            else:
                defaultname = ''
            templateDefName = 'template_'+tmplate+'.cif'
            if not CIFobj: # copying a template
                for pth in [os.getcwd()]+sys.path:
                    fil = os.path.join(pth,defaultname)
                    if os.path.exists(fil) and defaultname: break
                else:
                    for pth in sys.path:
                        fil = os.path.join(pth,templateDefName)
                        if os.path.exists(fil): break
                    else:
                        print(CIFobj+' not found in path!')
                        return
                fp = open(fil,'r')
                txt = fp.read()
                fp.close()
            elif type(CIFobj) is not list and type(CIFobj) is not tuple:
                if not os.path.exists(CIFobj):
                    print("Error: requested template file has disappeared: "+CIFobj)
                    return
                fp = open(CIFobj,'r')
                txt = fp.read()
                fp.close()
            else:
                txt = dict2CIF(CIFobj[0],CIFobj[1]).WriteOut()
            # remove the PyCifRW header, if present
            #if txt.find('PyCifRW') > -1 and txt.find('data_') > -1:
            txt = "# GSAS-II edited template follows "+txt[txt.index("data_")+5:]
            #txt = txt.replace('data_','#')
            WriteCIFitem(txt)

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
                        aname = str(phasedict['pId'])+':'+str(Histogram['hId'])+':'+item
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
                name = '%sBack;%d'%(hfx,i)
                sig = self.sigDict.get(name,-0.009)
                if len(l) > 60:
                    txt += l + '\n'
                    l = '    '
                l += G2mth.ValEsd(v,sig)+', '
            txt += l
            if bkgdict['nDebye']:
                txt += '\n  Background Debye function parameters: A, R, U:'
                names = ['A;','R;','U;']
                for i in range(bkgdict['nDebye']):
                    txt += '\n    '
                    for j in range(3):
                        name = hfx+'Debye'+names[j]+str(i)
                        sig = self.sigDict.get(name,-0.009)
                        txt += G2mth.ValEsd(bkgdict['debyeTerms'][i][2*j],sig)+', '
            if bkgdict['nPeaks']:
                txt += '\n  Background peak parameters: pos, int, sig, gam:'
                names = ['pos;','int;','sig;','gam;']
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
                s += '  peak variance(Gauss) = Utan(Th)^2^+Vtan(Th)+W:\n'
                s += '  peak HW(Lorentz) = X/cos(Th)+Ytan(Th); SH/L = S/L+H/L\n'
                s += '  U, V, W in (centideg)^2^, X & Y in centideg\n    '
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
                if len(self.powderDict) > 1:
                    if s:
                        s += '\n'
                    else:
                        s += '  Crystallite size model "%s" for %s (microns)\n  '%(size[0],phasenam)
                    s += '  Parameters for histogram #'+str(hId)+' '+str(histogram)+'\n'
                else:
                    s += '  Crystallite size model "%s" for %s (microns)\n  '%(size[0],phasenam)
                
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
                s += '\n  Mustrain model "%s" for %s (10^6^)\n  '%(mustrain[0],phasenam)
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
                s1 = '  \n  Macrostrain parameters: '
                names = G2spc.HStrainNames(SGData)
                for name in names:
                    s1 += name+', '
                s1 += '\n    '
                macrostrain = False
                for i in range(len(names)):
                    name = phfx+name[i]
                    sig = self.sigDict.get(name,-0.009)
                    s1 += G2mth.ValEsd(hstrain[0][i],sig)+', '
                    if hstrain[0][i]: macrostrain = True
                if macrostrain:
                    s += s1 + '\n'
                    # show revised lattice parameters here someday
                else:
                    s += '\n'
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
                         '\n   _atom_site_label'+
                         '\n   _atom_site_type_symbol'+
                         '\n   _atom_site_fract_x'+
                         '\n   _atom_site_fract_y'+
                         '\n   _atom_site_fract_z'+
                         '\n   _atom_site_occupancy'+
                         '\n   _atom_site_adp_type'+
                         '\n   _atom_site_U_iso_or_equiv'+
                         '\n   _atom_site_symmetry_multiplicity')

            varnames = {cx:'Ax',cx+1:'Ay',cx+2:'Az',cx+3:'Afrac',
                        cia+1:'AUiso',cia+2:'AU11',cia+3:'AU22',cia+4:'AU33',
                        cia+5:'AU12',cia+6:'AU13',cia+7:'AU23'}
            self.labellist = []
            
            pfx = str(phasedict['pId'])+'::'
            # loop over all atoms
            naniso = 0
            for i,at in enumerate(Atoms):
                if phasedict['General']['Type'] == 'macromolecular':
                    label = '%s_%s_%s_%s'%(at[ct-1],at[ct-3],at[ct-4],at[ct-2])
                    s = PutInCol(MakeUniqueLabel(label,self.labellist),15) # label
                else:
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
                for j in (cx,cx+1,cx+2,cx+3,cia,cia+1):
                    if j in (cx,cx+1,cx+2):
                        dig = 11
                        sigdig = -0.00009
                    else:
                        dig = 10
                        sigdig = -0.009
                    if j == cia:
                        s += adp
                    else:
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
                s += PutInCol(at[cs+1],3)
                WriteCIFitem(s)
            if naniso == 0: return
            # now loop over aniso atoms
            WriteCIFitem('\nloop_' + '\n   _atom_site_aniso_label' + 
                         '\n   _atom_site_aniso_U_11' + '\n   _atom_site_aniso_U_12' +
                         '\n   _atom_site_aniso_U_13' + '\n   _atom_site_aniso_U_22' +
                         '\n   _atom_site_aniso_U_23' + '\n   _atom_site_aniso_U_33')
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
                    print('Error: No mass found for atom type '+at[ct])
                    print('Will not compute cell contents for phase '+phasenam)
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
            # blocks  (for all atoms in all appropriate phases) - an example?:
#loop_
#    _atom_type_symbol
#    _atom_type_description
#    _atom_type_scat_dispersion_real
#    _atom_type_scat_dispersion_imag
#    _atom_type_scat_source
#    'C' 'C' 0.0033 0.0016
#                         'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
#    'H' 'H' 0.0000 0.0000
#                         'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
#    'P' 'P' 0.1023 0.0942
#                         'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
#    'Cl' 'Cl' 0.1484 0.1585
#                         'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
#    'Cu' 'Cu' 0.3201 1.2651
#                         'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'

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
            # create a dict for storing Pub flag for bonds/angles, if needed
            if phasedict['General'].get("DisAglHideFlag") is None:
                phasedict['General']["DisAglHideFlag"] = {}
            DisAngSel = phasedict['General']["DisAglHideFlag"]
            cx,ct,cs,cia = phasedict['General']['AtomPtrs']
            cn = ct-1
            fpfx = str(phasedict['pId'])+'::Afrac:'        
            cfrac = cx+3
            DisAglData = {}
            # create a list of atoms, but skip atoms with zero occupancy
            xyz = []
            fpfx = str(phasedict['pId'])+'::Afrac:'        
            for i,atom in enumerate(Atoms):
                if self.parmDict.get(fpfx+str(i),atom[cfrac]) == 0.0: continue
                xyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
            if 'DisAglCtls' not in generalData:
                # should not happen, since DisAglDialog should be called
                # for all phases before getting here
                dlg = G2gd.DisAglDialog(
                    self.G2frame,
                    {},
                    generalData)
                if dlg.ShowModal() == wx.ID_OK:
                    generalData['DisAglCtls'] = dlg.GetData()
                else:
                    dlg.Destroy()
                    return
                dlg.Destroy()
            DisAglData['OrigAtoms'] = xyz
            DisAglData['TargAtoms'] = xyz
            SymOpList,offsetList,symOpList,G2oprList = G2spc.AllOps(
                generalData['SGData'])

            xpandSGdata = generalData['SGData'].copy()
            xpandSGdata.update({'SGOps':symOpList,
                                'SGInv':False,
                                'SGLatt':'P',
                                'SGCen':np.array([[0, 0, 0]]),})
            DisAglData['SGData'] = xpandSGdata

            DisAglData['Cell'] = generalData['Cell'][1:] #+ volume
            if 'pId' in phasedict:
                DisAglData['pId'] = phasedict['pId']
                DisAglData['covData'] = self.OverallParms['Covariance']
            try:
                AtomLabels,DistArray,AngArray = G2stMn.RetDistAngle(
                    generalData['DisAglCtls'],
                    DisAglData)
            except KeyError:        # inside DistAngle for missing atom types in DisAglCtls
                print('**** ERROR - try again but do "Reset" to fill in missing atom types ****')
                    
            # loop over interatomic distances for this phase
            WriteCIFitem('\n# MOLECULAR GEOMETRY')
            WriteCIFitem('loop_' + 
                         '\n   _geom_bond_atom_site_label_1' +
                         '\n   _geom_bond_atom_site_label_2' + 
                         '\n   _geom_bond_distance' + 
                         '\n   _geom_bond_site_symmetry_1' + 
                         '\n   _geom_bond_site_symmetry_2' + 
                         '\n   _geom_bond_publ_flag')

            for i in sorted(AtomLabels.keys()):
                Dist = DistArray[i]
                for D in Dist:
                    line = '  '+PutInCol(AtomLabels[i],6)+PutInCol(AtomLabels[D[0]],6)
                    sig = D[4]
                    if sig == 0: sig = -0.00009
                    line += PutInCol(G2mth.ValEsd(D[3],sig,True),10)
                    line += "  1_555 "
                    line += " {:3d}_".format(D[2])
                    for d in D[1]:
                        line += "{:1d}".format(d+5)
                    if DisAngSel.get((i,tuple(D[0:3]))):
                        line += " no"
                    else:
                        line += " yes"
                    WriteCIFitem(line)

            # loop over interatomic angles for this phase
            WriteCIFitem('\nloop_' + 
                         '\n   _geom_angle_atom_site_label_1' + 
                         '\n   _geom_angle_atom_site_label_2' + 
                         '\n   _geom_angle_atom_site_label_3' + 
                         '\n   _geom_angle' + 
                         '\n   _geom_angle_site_symmetry_1' +
                         '\n   _geom_angle_site_symmetry_2' + 
                         '\n   _geom_angle_site_symmetry_3' + 
                         '\n   _geom_angle_publ_flag')

            for i in sorted(AtomLabels.keys()):
                Dist = DistArray[i]
                for k,j,tup in AngArray[i]:
                    Dj = Dist[j]
                    Dk = Dist[k]
                    line = '  '+PutInCol(AtomLabels[Dj[0]],6)+PutInCol(AtomLabels[i],6)+PutInCol(AtomLabels[Dk[0]],6)
                    sig = tup[1]
                    if sig == 0: sig = -0.009
                    line += PutInCol(G2mth.ValEsd(tup[0],sig,True),10)
                    line += " {:3d}_".format(Dj[2])
                    for d in Dj[1]:
                        line += "{:1d}".format(d+5)
                    line += "  1_555 "
                    line += " {:3d}_".format(Dk[2])
                    for d in Dk[1]:
                        line += "{:1d}".format(d+5)
                    key = (tuple(Dk[0:3]),i,tuple(Dj[0:3]))
                    if DisAngSel.get(key):
                        line += " no"
                    else:
                        line += " yes"
                    WriteCIFitem(line)

        def WritePhaseInfo(phasenam):
            'Write out the phase information for the selected phase'
            WriteCIFitem('\n# phase info for '+str(phasenam) + ' follows')
            phasedict = self.Phases[phasenam] # pointer to current phase info            
            WriteCIFitem('_pd_phase_name', phasenam)
            cellList,cellSig = self.GetCell(phasenam)
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
            WriteCIFitem('loop_\n    _space_group_symop_id\n    _space_group_symop_operation_xyz')
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
                            print("Internal error: no block for data. Phase "+str(
                                phasenam)+" histogram "+str(hist))
                            histlist = []
                            break
                        histlist.append(blockid)

                if len(histlist) == 0:
                    WriteCIFitem('# Note: phase has no associated data')

            # report atom params
            if phasedict['General']['Type'] in ['nuclear','macromolecular']:        #this needs macromolecular variant, etc!
                WriteAtomsNuclear(phasenam)
            else:
                raise Exception,"no export for "+str(phasedict['General']['Type'])+" coordinates implemented"
            # report cell contents
            WriteComposition(phasenam)
            if not self.quickmode and phasedict['General']['Type'] == 'nuclear':      # report distances and angles
                WriteDistances(phasenam,SymOpList,offsetList,symOpList,G2oprList)
            if 'Map' in phasedict['General'] and 'minmax' in phasedict['General']['Map']:
                WriteCIFitem('\n# Difference density results')
                MinMax = phasedict['General']['Map']['minmax']
                WriteCIFitem('_refine_diff_density_max',G2mth.ValEsd(MinMax[1],-0.009))
                WriteCIFitem('_refine_diff_density_min',G2mth.ValEsd(MinMax[0],-0.009))
                
        def Yfmt(ndec,val):
            'Format intensity values'
            out = ("{:."+str(ndec)+"f}").format(val)
            out = out.rstrip('0')  # strip zeros to right of decimal
            return out.rstrip('.')  # and decimal place when not needed
            
        def WriteReflStat(refcount,hklmin,hklmax,dmin,dmax,nRefSets=1):
            'Write reflection statistics'
            WriteCIFitem('_reflns_number_total', str(refcount))
            if hklmin is not None and nRefSets == 1: # hkl range has no meaning with multiple phases
                WriteCIFitem('_reflns_limit_h_min', str(int(hklmin[0])))
                WriteCIFitem('_reflns_limit_h_max', str(int(hklmax[0])))
                WriteCIFitem('_reflns_limit_k_min', str(int(hklmin[1])))
                WriteCIFitem('_reflns_limit_k_max', str(int(hklmax[1])))
                WriteCIFitem('_reflns_limit_l_min', str(int(hklmin[2])))
                WriteCIFitem('_reflns_limit_l_max', str(int(hklmax[2])))
            if hklmin is not None:
                WriteCIFitem('_reflns_d_resolution_low  ', G2mth.ValEsd(dmax,-0.009))
                WriteCIFitem('_reflns_d_resolution_high ', G2mth.ValEsd(dmin,-0.009))

        def WritePowderData(histlbl):
            'Write out the selected powder diffraction histogram info'
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
                WriteCIFitem('_diffrn_radiation_type','K\\a~1,2~')
                WriteCIFitem('loop_' + 
                             '\n   _diffrn_radiation_wavelength' +
                             '\n   _diffrn_radiation_wavelength_wt' + 
                             '\n   _diffrn_radiation_wavelength_id')
                WriteCIFitem('  ' + PutInCol(G2mth.ValEsd(lam1,slam1),15)+
                             PutInCol('1.0',15) + 
                             PutInCol('1',5))
                WriteCIFitem('  ' + PutInCol(G2mth.ValEsd(lam2,slam2),15)+
                             PutInCol(G2mth.ValEsd(ratio,sratio),15)+
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
                                 '\n   _pd_phase_id' + 
                                 '\n   _pd_phase_block_id' + 
                                 '\n   _pd_phase_mass_%')
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
                    WriteCIFitem('loop_' +
                                 '\n   _gsas_proc_phase_R_F_factor' +
                                 '\n   _gsas_proc_phase_R_Fsqd_factor' +
                                 '\n   _gsas_proc_phase_id' +
                                 '\n   _gsas_proc_phase_block_id')
                    for phasenam in phasebyhistDict.get(histlbl):
                        pfx = str(self.Phases[phasenam]['pId'])+':'+str(hId)+':'
                        WriteCIFitem(
                            '  '+
                            '  '+G2mth.ValEsd(histblk[pfx+'Rf']/100.,-.00009) +
                            '  '+G2mth.ValEsd(histblk[pfx+'Rf^2']/100.,-.00009)+
                            '  '+str(self.Phases[phasenam]['pId'])+
                            '  '+datablockidDict[phasenam]
                            )
            else:
                # single phase in this histogram
                pfx = '0:'+str(hId)+':'
                WriteCIFitem('_refine_ls_R_F_factor      ','%.5f'%(histblk[pfx+'Rf']/100.))
                WriteCIFitem('_refine_ls_R_Fsqd_factor   ','%.5f'%(histblk[pfx+'Rf^2']/100.))
                
            WriteCIFitem('_pd_proc_ls_prof_R_factor   ','%.5f'%(histblk['R']/100.))
            WriteCIFitem('_pd_proc_ls_prof_wR_factor  ','%.5f'%(histblk['wR']/100.))
            WriteCIFitem('_gsas_proc_ls_prof_R_B_factor ','%.5f'%(histblk['Rb']/100.))
            WriteCIFitem('_gsas_proc_ls_prof_wR_B_factor','%.5f'%(histblk['wRb']/100.))
            WriteCIFitem('_pd_proc_ls_prof_wR_expected','%.5f'%(histblk['wRmin']/100.))

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
            # data collection parameters for the powder dataset
            
            temperature = histblk['Sample Parameters'].get('Temperature') # G2 uses K
            if not temperature:
                T = '?'
            else:
                T = G2mth.ValEsd(temperature,-0.009,True) # CIF uses K
            WriteCIFitem('_diffrn_ambient_temperature',T)

            pressure = histblk['Sample Parameters'].get('Pressure') #G2 uses mega-Pascal
            if not pressure:
                P = '?'
            else:
                P = G2mth.ValEsd(pressure*1000,-0.09,True) # CIF uses kilopascal
            WriteCIFitem('_diffrn_ambient_pressure',P)

            WriteCIFitem('\n# STRUCTURE FACTOR TABLE')            
            # compute maximum intensity reflection
            Imax = 0
            for phasenam in histblk['Reflection Lists']:
                scale = self.Phases[phasenam]['Histograms'][histlbl]['Scale'][0]
                refList = np.asarray(histblk['Reflection Lists'][phasenam]['RefList'])
                I100 = scale*refList.T[8]*refList.T[11]
                #Icorr = np.array([refl[13] for refl in histblk['Reflection Lists'][phasenam]])[0]
                #FO2 = np.array([refl[8] for refl in histblk['Reflection Lists'][phasenam]])
                #I100 = scale*FO2*Icorr
                Imax = max(Imax,max(I100))

            WriteCIFitem('loop_')
            if len(histblk['Reflection Lists'].keys()) > 1:
                WriteCIFitem('   _pd_refln_phase_id')
            WriteCIFitem('   ' + refprx + 'index_h' + 
                         '\n   ' + refprx + 'index_k' + 
                         '\n   ' + refprx + 'index_l' + 
                         '\n   ' + refprx + 'F_squared_meas' + 
                         '\n   ' + refprx + 'F_squared_calc' + 
                         '\n   ' + refprx + 'phase_calc' + 
                         '\n   _pd_refln_d_spacing')
            if Imax > 0:
                WriteCIFitem('   _gsas_i100_meas')

            refcount = 0
            hklmin = None
            hklmax = None
            dmax = None
            dmin = None
            for phasenam in histblk['Reflection Lists']:
                scale = self.Phases[phasenam]['Histograms'][histlbl]['Scale'][0]
                phaseid = self.Phases[phasenam]['pId']
                refcount += len(histblk['Reflection Lists'][phasenam]['RefList'])
                for j,ref in enumerate(histblk['Reflection Lists'][phasenam]['RefList']):
                    if DEBUG:
                        print('DEBUG: skipping reflection list')
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
                        s += PutInCol(G2mth.ValEsd(100.*I100[j]/Imax,-0.09),6)
                    WriteCIFitem("  "+s)

            WriteReflStat(refcount,hklmin,hklmax,dmin,dmax,len(histblk['Reflection Lists']))
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
            #            WriteCIFitem('   _pd_proc_d_spacing') # need easy way to get this
            if not fixedstep:
                if zero:
                    WriteCIFitem('   _pd_proc_2theta_corrected')
                else:
                    WriteCIFitem('   _pd_meas_2theta_scan')
            # at least for now, always report weights.
            #if countsdata:
            #    WriteCIFitem('   _pd_meas_counts_total')
            #else:
            WriteCIFitem('   _pd_meas_intensity_total')
            WriteCIFitem('   _pd_calc_intensity_total')
            WriteCIFitem('   _pd_proc_intensity_bkg_calc')
            WriteCIFitem('   _pd_proc_ls_weight')
            maxY = max(histblk['Data'][1].max(),histblk['Data'][3].max())
            if maxY < 0: maxY *= -10 # this should never happen, but...
            ndec = max(0,10-int(np.log10(maxY))-1) # 10 sig figs should be enough
            maxSU = histblk['Data'][2].max()
            if maxSU < 0: maxSU *= -1 # this should never happen, but...
            ndecSU = max(0,8-int(np.log10(maxSU))-1) # 8 sig figs should be enough
            lowlim,highlim = histblk['Limits'][1]

            if DEBUG:
                print('DEBUG: skipping profile list')
            else:    
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
            'Write out the selected single crystal histogram info'
            histblk = self.Histograms[histlbl]

            #refprx = '_refln.' # mm
            refprx = '_refln_' # normal

            WriteCIFitem('\n# STRUCTURE FACTOR TABLE')            
            WriteCIFitem('loop_' + 
                         '\n   ' + refprx + 'index_h' + 
                         '\n   ' + refprx + 'index_k' + 
                         '\n   ' + refprx + 'index_l' +
                         '\n   ' + refprx + 'F_squared_meas' + 
                         '\n   ' + refprx + 'F_squared_sigma' + 
                         '\n   ' + refprx + 'F_squared_calc' + 
                         '\n   ' + refprx + 'phase_calc'
                         )

            hklmin = None
            hklmax = None
            dmax = None
            dmin = None
            refcount = len(histblk['Data']['RefList'])
            for ref in histblk['Data']['RefList']:
                s = "  "
                if hklmin is None:
                    hklmin = ref[0:3]
                    hklmax = ref[0:3]
                    dmax = dmin = ref[4]
                for i,hkl in enumerate(ref[0:3]):
                    hklmax[i] = max(hkl,hklmax[i])
                    hklmin[i] = min(hkl,hklmin[i])
                    s += PutInCol(int(hkl),4)
                import sys
                if ref[5] == 0.0:
                    s += PutInCol(G2mth.ValEsd(ref[8],0),12)
                    s += PutInCol('.',10)
                    s += PutInCol(G2mth.ValEsd(ref[9],0),12)
                    s += PutInCol(G2mth.ValEsd(ref[10],-0.9),7)
                else:
                    sig = ref[6] * ref[8] / ref[5]
                    s += PutInCol(G2mth.ValEsd(ref[8],-abs(sig/10)),12)
                    s += PutInCol(G2mth.ValEsd(sig,-abs(sig)/10.),10)
                    s += PutInCol(G2mth.ValEsd(ref[9],-abs(sig/10)),12)
                s += PutInCol(G2mth.ValEsd(ref[10],-0.9),7)
                dmax = max(dmax,ref[4])
                dmin = min(dmin,ref[4])
                WriteCIFitem(s)
            if not self.quickmode: # statistics only in a full CIF
                WriteReflStat(refcount,hklmin,hklmax,dmin,dmax)
                hId = histblk['hId']
                hfx = '0:'+str(hId)+':'
                phfx = '%d:%d:'%(0,hId)
                extType,extModel,extParms = self.Phases[phasenam]['Histograms'][histlbl]['Extinction']
                if extModel != 'None':
                    WriteCIFitem('# Extinction scaled by 1.e5')
                    WriteCIFitem('_refine_ls_extinction_method','Becker-Coppens %s %s'%(extModel,extType))
                    sig = -1.e-3
                    if extModel == 'Primary':
                        parm = extParms['Ep'][0]*1.e5
                        if extParms['Ep'][1]:
                            sig = self.sigDict[phfx+'Ep']*1.e5
                        text = G2mth.ValEsd(parm,sig)
                    elif extModel == 'Secondary Type I':
                        parm = extParms['Eg'][0]*1.e5
                        if extParms['Eg'][1]:
                            sig = self.sigDict[phfx+'Eg']*1.e5
                        text = G2mth.ValEsd(parm,sig)
                    elif extModel == 'Secondary Type II':
                        parm = extParms['Es'][0]*1.e5
                        if extParms['Es'][1]:
                            sig = self.sigDict[phfx+'Es']*1.e5
                        text = G2mth.ValEsd(parm,sig)
                    elif extModel == 'Secondary Type I & II':
                        parm = extParms['Eg'][0]*1.e5
                        if extParms['Es'][1]:
                            sig = self.sigDict[phfx+'Es']*1.e5
                        text = G2mth.ValEsd(parm,sig)
                        sig = -1.0e-3
                        parm = extParms['Es'][0]*1.e5
                        if extParms['Es'][1]:
                            sig = self.sigDict[phfx+'Es']*1.e5
                        text += G2mth.ValEsd(parm,sig)
                    WriteCIFitem('_refine_ls_extinction_coef',text)
                    WriteCIFitem('_refine_ls_extinction_expression','Becker & Coppens (1974). Acta Cryst. A30, 129-147')

                WriteCIFitem('_refine_ls_wR_factor_gt    ','%.4f'%(histblk['wR']/100.))
                WriteCIFitem('_refine_ls_R_factor_gt     ','%.4f'%(histblk[hfx+'Rf']/100.))
                WriteCIFitem('_refine_ls_R_Fsqd_factor   ','%.4f'%(histblk[hfx+'Rf^2']/100.))
        def EditAuthor(event=None):
            'dialog to edit the CIF author info'
            'Edit the CIF author name'
            dlg = G2gd.SingleStringDialog(self.G2frame,
                                          'Get CIF Author',
                                          'Provide CIF Author name (Last, First)',
                                          value=self.author)
            if not dlg.Show():
                dlg.Destroy()
                return False  # cancel was pressed
            self.author = dlg.GetValue()
            dlg.Destroy()
            try:
                self.OverallParms['Controls']["Author"] = self.author # save for future
            except KeyError:
                pass
            return True
        def EditInstNames(event=None):
            'Provide a dialog for editing instrument names'
            dictlist = []
            keylist = []
            lbllist = []
            for hist in self.Histograms:
                if hist.startswith("PWDR"): 
                    key2 = "Sample Parameters"
                    d = self.Histograms[hist][key2]
                elif hist.startswith("HKLF"): 
                    key2 = "Instrument Parameters"
                    d = self.Histograms[hist][key2][0]
                    
                lbllist.append(hist)
                dictlist.append(d)
                keylist.append('InstrName')
                instrname = d.get('InstrName')
                if instrname is None:
                    d['InstrName'] = ''
            return G2ctrls.CallScrolledMultiEditor(
                self.G2frame,dictlist,keylist,
                prelbl=range(1,len(dictlist)+1),
                postlbl=lbllist,
                title='Instrument names',
                header="Edit instrument names. Note that a non-blank\nname is required for all histograms",
                CopyButton=True)
            
        def EditRanges(event):
            '''Edit the bond distance/angle search range; phase is determined from
            a pointer placed in the button object (.phasedict) that references the
            phase dictionary
            '''
            but = event.GetEventObject()
            phasedict = but.phasedict
            dlg = G2gd.DisAglDialog(
                self.G2frame,
                phasedict['General']['DisAglCtls'], # edited 
                phasedict['General'], # defaults
                )
            if dlg.ShowModal() == wx.ID_OK:
                phasedict['General']['DisAglCtls'] = dlg.GetData()
            dlg.Destroy()
            
        def EditCIFDefaults():
            '''Fills the CIF Defaults window with controls for editing various CIF export
            parameters (mostly related to templates).
            '''
            self.cifdefs.DestroyChildren()
            self.cifdefs.SetTitle('Edit CIF settings')
            vbox = wx.BoxSizer(wx.VERTICAL)
            vbox.Add(wx.StaticText(self.cifdefs, wx.ID_ANY,'Creating file '+str(self.filename)))
            but = wx.Button(self.cifdefs, wx.ID_ANY,'Edit CIF Author')
            but.Bind(wx.EVT_BUTTON,EditAuthor)
            vbox.Add(but,0,wx.ALIGN_CENTER,3)
            but = wx.Button(self.cifdefs, wx.ID_ANY,'Edit Instrument Name(s)')
            but.Bind(wx.EVT_BUTTON,EditInstNames)
            vbox.Add(but,0,wx.ALIGN_CENTER,3)
            cpnl = wxscroll.ScrolledPanel(self.cifdefs,size=(300,300))
            cbox = wx.BoxSizer(wx.VERTICAL)
            G2gd.HorizontalLine(cbox,cpnl)          
            cbox.Add(
                CIFtemplateSelect(self.cifdefs,
                                  cpnl,'publ',self.OverallParms['Controls'],
                                  EditCIFDefaults,
                                  "Publication (overall) template",
                                  ),
                0,wx.EXPAND|wx.ALIGN_LEFT|wx.ALL)
            for phasenam in sorted(self.Phases.keys()):
                G2gd.HorizontalLine(cbox,cpnl)          
                title = 'Phase '+phasenam
                phasedict = self.Phases[phasenam] # pointer to current phase info            
                cbox.Add(
                    CIFtemplateSelect(self.cifdefs,
                                      cpnl,'phase',phasedict['General'],
                                      EditCIFDefaults,
                                      title,
                                      phasenam),
                    0,wx.EXPAND|wx.ALIGN_LEFT|wx.ALL)
                cpnl.SetSizer(cbox)
                if phasedict['General']['Type'] == 'nuclear':  
                    but = wx.Button(cpnl, wx.ID_ANY,'Edit distance/angle ranges')
                    cbox.Add(but,0,wx.ALIGN_LEFT,0)
                    cbox.Add((-1,2))
                    but.phasedict = self.Phases[phasenam]  # set a pointer to current phase info
                    but.Bind(wx.EVT_BUTTON,EditRanges)     # phase bond/angle ranges
                    but = wx.Button(cpnl, wx.ID_ANY,'Set distance/angle publication flags')
                    but.phase = phasenam  # set a pointer to current phase info     
                    but.Bind(wx.EVT_BUTTON,SelectDisAglFlags)     # phase bond/angle ranges
                    cbox.Add(but,0,wx.ALIGN_LEFT,0)
                cbox.Add((-1,2))
            for i in sorted(self.powderDict.keys()):
                G2gd.HorizontalLine(cbox,cpnl)          
                hist = self.powderDict[i]
                histblk = self.Histograms[hist]
                title = 'Powder dataset '+hist[5:]
                cbox.Add(
                    CIFtemplateSelect(self.cifdefs,
                                      cpnl,'powder',histblk["Sample Parameters"],
                                      EditCIFDefaults,
                                      title,
                                      histblk["Sample Parameters"]['InstrName']),
                    0,wx.EXPAND|wx.ALIGN_LEFT|wx.ALL)
            for i in sorted(self.xtalDict.keys()):
                G2gd.HorizontalLine(cbox,cpnl)          
                hist = self.xtalDict[i]
                histblk = self.Histograms[hist]
                title = 'Single Xtal dataset '+hist[5:]
                cbox.Add(
                    CIFtemplateSelect(self.cifdefs,
                                      cpnl,'single',histblk["Instrument Parameters"][0],
                                      EditCIFDefaults,
                                      title,
                                      histblk["Instrument Parameters"][0]['InstrName']),
                    0,wx.EXPAND|wx.ALIGN_LEFT|wx.ALL)
            cpnl.SetSizer(cbox)
            cpnl.SetAutoLayout(1)
            cpnl.SetupScrolling()
            #cpnl.Bind(rw.EVT_RW_LAYOUT_NEEDED, self.OnLayoutNeeded) # needed if sizes change
            cpnl.Layout()

            vbox.Add(cpnl, 1, wx.ALIGN_LEFT|wx.ALL|wx.EXPAND, 0)
            btnsizer = wx.StdDialogButtonSizer()
            btn = wx.Button(self.cifdefs, wx.ID_OK, "Create CIF")
            btn.SetDefault()
            btnsizer.AddButton(btn)
            btn = wx.Button(self.cifdefs, wx.ID_CANCEL)
            btnsizer.AddButton(btn)
            btnsizer.Realize()
            vbox.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
            self.cifdefs.SetSizer(vbox)
            vbox.Fit(self.cifdefs)
            self.cifdefs.Layout()
            
        def OnToggleButton(event):
            'Respond to press of ToggleButton in SelectDisAglFlags'
            but = event.GetEventObject()
            if but.GetValue():                    
                but.DisAglSel[but.key] = True
            else:
                try:
                    del but.DisAglSel[but.key]
                except KeyError:
                    pass
        def keepTrue(event):
            event.GetEventObject().SetValue(True)
        def keepFalse(event):
            event.GetEventObject().SetValue(False)
                
        def SelectDisAglFlags(event):
            'Select Distance/Angle use flags for the selected phase'
            phasenam = event.GetEventObject().phase
            phasedict = self.Phases[phasenam]
            SymOpList,offsetList,symOpList,G2oprList = G2spc.AllOps(phasedict['General']['SGData'])
            generalData = phasedict['General']
            # create a dict for storing Pub flag for bonds/angles, if needed
            if phasedict['General'].get("DisAglHideFlag") is None:
                phasedict['General']["DisAglHideFlag"] = {}
            DisAngSel = phasedict['General']["DisAglHideFlag"]

            cx,ct,cs,cia = phasedict['General']['AtomPtrs']
            cn = ct-1
            cfrac = cx+3
            DisAglData = {}
            # create a list of atoms, but skip atoms with zero occupancy
            xyz = []
            fpfx = str(phasedict['pId'])+'::Afrac:'        
            for i,atom in enumerate(phasedict['Atoms']):
                if self.parmDict.get(fpfx+str(i),atom[cfrac]) == 0.0: continue
                xyz.append([i,]+atom[cn:cn+2]+atom[cx:cx+3])
            if 'DisAglCtls' not in generalData:
                # should not be used, since DisAglDialog should be called
                # for all phases before getting here
                dlg = G2gd.DisAglDialog(
                    self.cifdefs,
                    {},
                    generalData)
                if dlg.ShowModal() == wx.ID_OK:
                    generalData['DisAglCtls'] = dlg.GetData()
                else:
                    dlg.Destroy()
                    return
                dlg.Destroy()
            dlg = wx.Dialog(
                self.G2frame,
                style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
            vbox = wx.BoxSizer(wx.VERTICAL)
            txt = wx.StaticText(dlg,wx.ID_ANY,'Searching distances for phase '+phasenam
                                +'\nPlease wait...')
            vbox.Add(txt,0,wx.ALL|wx.EXPAND)
            dlg.SetSizer(vbox)
            dlg.CenterOnParent()
            dlg.Show() # post "please wait"
            wx.BeginBusyCursor() # and change cursor

            DisAglData['OrigAtoms'] = xyz
            DisAglData['TargAtoms'] = xyz
            SymOpList,offsetList,symOpList,G2oprList = G2spc.AllOps(
                generalData['SGData'])

            xpandSGdata = generalData['SGData'].copy()
            xpandSGdata.update({'SGOps':symOpList,
                                'SGInv':False,
                                'SGLatt':'P',
                                'SGCen':np.array([[0, 0, 0]]),})
            DisAglData['SGData'] = xpandSGdata

            DisAglData['Cell'] = generalData['Cell'][1:] #+ volume
            if 'pId' in phasedict:
                DisAglData['pId'] = phasedict['pId']
                DisAglData['covData'] = self.OverallParms['Covariance']
            try:
                AtomLabels,DistArray,AngArray = G2stMn.RetDistAngle(
                    generalData['DisAglCtls'],
                    DisAglData)
            except KeyError:        # inside DistAngle for missing atom types in DisAglCtls
                print('**** ERROR - try again but do "Reset" to fill in missing atom types ****')
            wx.EndBusyCursor()
            txt.SetLabel('Set publication flags for distances and angles in\nphase '+phasenam)
            vbox.Add((5,5)) 
            vbox.Add(wx.StaticText(dlg,wx.ID_ANY,
                                   'The default is to flag all distances and angles as to be'+
                                   '\npublished. Change this by pressing appropriate buttons.'),
                     0,wx.ALL|wx.EXPAND)
            hbox = wx.BoxSizer(wx.HORIZONTAL)
            vbox.Add(hbox)
            hbox.Add(wx.StaticText(dlg,wx.ID_ANY,'Button appearance: '))
            but = wx.ToggleButton(dlg,wx.ID_ANY,'Publish')
            but.Bind(wx.EVT_TOGGLEBUTTON,keepFalse)
            hbox.Add(but)
            but = wx.ToggleButton(dlg,wx.ID_ANY,"Don't publish")
            but.Bind(wx.EVT_TOGGLEBUTTON,keepTrue)
            hbox.Add(but)
            but.SetValue(True)
            G2gd.HorizontalLine(vbox,dlg)          
            
            cpnl = wxscroll.ScrolledPanel(dlg,size=(400,300))
            cbox = wx.BoxSizer(wx.VERTICAL)
            for c in sorted(DistArray):
                karr = []
                UsedCols = {}
                cbox.Add(wx.StaticText(cpnl,wx.ID_ANY,
                                   'distances to/angles around atom '+AtomLabels[c]))
                #dbox = wx.GridBagSizer(hgap=5)
                dbox = wx.GridBagSizer()
                for i,D in enumerate(DistArray[c]):
                    karr.append(tuple(D[0:3]))
                    val = "{:.2f}".format(D[3])
                    sym = " [{:d} {:d} {:d}]".format(*D[1]) + " #{:d}".format(D[2])
                    dbox.Add(wx.StaticText(cpnl,wx.ID_ANY,AtomLabels[D[0]]),
                             (i+1,0)
                             )                    
                    dbox.Add(wx.StaticText(cpnl,wx.ID_ANY,sym),
                             (i+1,1)
                             )
                    but = wx.ToggleButton(cpnl,wx.ID_ANY,val)
                    but.key = (c,karr[-1])
                    but.DisAglSel = DisAngSel
                    if DisAngSel.get(but.key): but.SetValue(True)
                    but.Bind(wx.EVT_TOGGLEBUTTON,OnToggleButton)
                    dbox.Add(but,(i+1,2),border=1)
                for i,D in enumerate(AngArray[c]):
                    val = "{:.1f}".format(D[2][0])
                    but = wx.ToggleButton(cpnl,wx.ID_ANY,val)
                    but.key = (karr[D[0]],c,karr[D[1]])
                    but.DisAglSel = DisAngSel
                    if DisAngSel.get(but.key): but.SetValue(True)
                    but.Bind(wx.EVT_TOGGLEBUTTON,OnToggleButton)
                    dbox.Add(but,(D[0]+1,D[1]+3),border=1)
                    UsedCols[D[1]+3] = True
                for i,D in enumerate(DistArray[c][:-1]): # label columns that are used
                    if UsedCols.get(i+3):
                        dbox.Add(wx.StaticText(cpnl,wx.ID_ANY,AtomLabels[D[0]]),
                                 (0,i+3),
                                 flag=wx.ALIGN_CENTER
                                 )
                dbox.Add(wx.StaticText(cpnl,wx.ID_ANY,'distance'),
                                 (0,2),
                                 flag=wx.ALIGN_CENTER
                                 )
                cbox.Add(dbox)
                G2gd.HorizontalLine(cbox,cpnl)          
            cpnl.SetSizer(cbox)
            cpnl.SetAutoLayout(1)
            cpnl.SetupScrolling()
            #cpnl.Bind(rw.EVT_RW_LAYOUT_NEEDED, self.OnLayoutNeeded) # needed if sizes change
            cpnl.Layout()

            vbox.Add(cpnl, 1, wx.ALIGN_LEFT|wx.ALL|wx.EXPAND, 0)

            btnsizer = wx.StdDialogButtonSizer()
            btn = wx.Button(dlg, wx.ID_OK, "Done")
            btn.SetDefault()
            btnsizer.AddButton(btn)
            btnsizer.Realize()
            vbox.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
            dlg.SetSizer(vbox)
            vbox.Fit(dlg)
            dlg.Layout()
            
            dlg.CenterOnParent()
            dlg.ShowModal()
            
#***** end of functions for export method =======================================
#=================================================================================

        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()
        if self.mode=='simple':
            if self.ExportSelect('ask'): return
        else:
            if self.ExportSelect('default'): return
        # Someday: get restraint & constraint info
        #restraintDict = self.OverallParms.get('Restraints',{})
        #for i in  self.OverallParms['Constraints']:
        #    print i
        #    for j in self.OverallParms['Constraints'][i]:
        #        print j

        self.CIFdate = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%dT%H:%M")
        # is there anything to export?
        if len(self.Phases) == len(self.powderDict) == len(self.xtalDict) == 0:
           self.G2frame.ErrorDialog(
               'Empty project',
               'Project does not contain any data or phases. Are they interconnected?')
           return
        # get the project file name
        self.CIFname = os.path.splitext(
            os.path.split(self.G2frame.GSASprojectfile)[1]
            )[0]
        self.CIFname = self.CIFname.replace(' ','')
        if not self.CIFname: # none defined & needed, save as GPX to get one
            self.G2frame.OnFileSaveas(None)
            if not self.G2frame.GSASprojectfile: return
            self.CIFname = os.path.splitext(
                os.path.split(self.G2frame.GSASprojectfile)[1]
                )[0]
            self.CIFname = self.CIFname.replace(' ','')
        # test for quick CIF mode or no data
        self.quickmode = False
        phasenam = None # include all phases
        if self.mode == "simple" and self.currentExportType == 'phase':
            if len(self.Phases) == 0: # this check is probably not needed
                self.G2frame.ErrorDialog(
                    'No phase present',
                    'Cannot create a coordinates CIF with no phases')
                return
            self.quickmode = True
            oneblock = True
            if len(self.Phases) > 1: # quick mode: get selected phase
                phasenam = self.phasenam[0]
        elif self.mode == "simple": # powder/single xtal data export
            self.quickmode = True
            oneblock = True
        # Project export: will this require a multiblock CIF?
        elif len(self.Phases) > 1:
            oneblock = False
        elif len(self.powderDict) + len(self.xtalDict) > 1:
            oneblock = False
        else: # one phase, one dataset, Full CIF
            oneblock = True

        # make sure required information is present
        # get CIF author name -- required for full CIFs
        try:
            self.author = self.OverallParms['Controls'].get("Author",'').strip()
        except KeyError:
            pass
        while not (self.author or self.quickmode):
            if not EditAuthor(): return
        self.shortauthorname = self.author.replace(',','').replace(' ','')[:20]

        # check there is an instrument name for every histogram
        self.ifPWDR = False
        self.ifHKLF = False
        if not self.quickmode:
            invalid = 0
            key3 = 'InstrName'
            for hist in self.Histograms:
                if hist.startswith("PWDR"):
                    self.ifPWDR = True
                    key2 = "Sample Parameters"
                    d = self.Histograms[hist][key2]
                elif hist.startswith("HKLF"):
                    self.ifHKLF = True
                    key2 = "Instrument Parameters"
                    d = self.Histograms[hist][key2][0]                    
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
                if not EditInstNames(): return
        if not self.quickmode:
            # check for a distance-angle range search range for each phase
            for phasenam in sorted(self.Phases.keys()):
                #i = self.Phases[phasenam]['pId']
                phasedict = self.Phases[phasenam] # pointer to current phase info            
                if 'DisAglCtls' not in phasedict['General']:
                    dlg = G2gd.DisAglDialog(
                        self.G2frame,
                        {},
                        phasedict['General'])
                    if dlg.ShowModal() == wx.ID_OK:
                        phasedict['General']['DisAglCtls'] = dlg.GetData()
                    else:
                        dlg.Destroy()
                        return
                    dlg.Destroy()
        if not self.quickmode:
            # check if temperature values & pressure are defaulted
            default = 0
            for hist in self.Histograms:
                if hist.startswith("PWDR"): 
                    key2 = "Sample Parameters"
                    T = self.Histograms[hist][key2].get('Temperature')
                    if not T:
                        default += 1
                    elif T == 300:
                        default += 1
                    P = self.Histograms[hist][key2].get('Pressure')
                    if not P:
                        default += 1
                    elif P == 1:
                        default += 1
            if default > 0:
                dlg = wx.MessageDialog(
                    self.G2frame,
                    'Temperature/Pressure values appear to be defaulted for some powder histograms (See Sample Parameters for each PWDR tree entry). Do you want to use those values?',
                    'Check T and P values',
                    wx.OK|wx.CANCEL)
                ret = dlg.ShowModal()
                dlg.Destroy()
                if ret != wx.ID_OK: return
        if oneblock and not self.quickmode:
            # select a dataset to use (there should only be one set in one block,
            # but take whatever comes 1st)
            for hist in self.Histograms:
                histblk = self.Histograms[hist]
                if hist.startswith("PWDR"): 
                    instnam = histblk["Sample Parameters"]['InstrName']
                    break # ignore all but 1st data histogram
                elif hist.startswith("HKLF"): 
                    instnam = histblk["Instrument Parameters"][0]['InstrName']
                    break # ignore all but 1st data histogram
        if not self.quickmode: # give the user a chance to edit all defaults
            self.cifdefs = wx.Dialog(
                self.G2frame,
                style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
            EditCIFDefaults()
            self.cifdefs.CenterOnParent()
            val = self.cifdefs.ShowModal()
            self.cifdefs.Destroy()
            if val != wx.ID_OK:
                return
        #======================================================================
        # Start writing the CIF - single block
        #======================================================================
        print('Writing CIF output to file '+str(self.filename)+"...")
        self.OpenFile()
        if self.currentExportType == 'single' or self.currentExportType == 'powder':
            #======Data only CIF (powder/xtal) ====================================
            hist = self.histnam[0]
            self.CIFname = hist[5:40].replace(' ','')
            WriteCIFitem('data_'+self.CIFname)
            if hist.startswith("PWDR"):
                WritePowderData(hist)
            elif hist.startswith("HKLF"):
                WriteSingleXtalData(hist)
            else:
                print "should not happen"
        elif self.quickmode:
            #====Phase only CIF ====================================================
            WriteCIFitem('data_'+self.CIFname)
            if phasenam is None: # if not already selected, select the first phase (should be one) 
                phasenam = self.Phases.keys()[0]
            #print 'phasenam',phasenam
            phaseblk = self.Phases[phasenam] # pointer to current phase info
            # report the phase info
            WritePhaseInfo(phasenam)
        elif oneblock:
            #====Single block, data & phase CIF ===================================
            WriteCIFitem('data_'+self.CIFname)
            if phasenam is None: # if not already selected, select the first phase (should be one) 
                phasenam = self.Phases.keys()[0]
            #print 'phasenam',phasenam
            phaseblk = self.Phases[phasenam] # pointer to current phase info
            instnam = instnam.replace(' ','')
            WriteCIFitem('_pd_block_id',
                         str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                         str(self.shortauthorname) + "|" + instnam)
            WriteAudit()
            writeCIFtemplate(self.OverallParms['Controls'],'publ') # overall (publication) template
            WriteOverall()
            writeCIFtemplate(self.Phases[phasenam]['General'],'phase',phasenam) # write phase template
            # report the phase info
            WritePhaseInfo(phasenam)
            if hist.startswith("PWDR"):
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
                histblk = self.Histograms[hist]["Sample Parameters"]
                writeCIFtemplate(histblk,'powder',histblk['InstrName']) # write powder template
                WritePowderData(hist)
            elif hist.startswith("HKLF"):
                histprm = self.Histograms[hist]["Instrument Parameters"][0]
                writeCIFtemplate(histprm,'single',histprm['InstrName']) # single crystal template
                WriteSingleXtalData(hist)
        else:
            #=== multiblock: multiple phases and/or histograms ====================
            nsteps = 1 + len(self.Phases) + len(self.powderDict) + len(self.xtalDict)
            try:
                dlg = wx.ProgressDialog('CIF progress','starting',nsteps,parent=self.G2frame)
                Size = dlg.GetSize()
                Size = (int(Size[0]*3),Size[1]) # increase size along x
                dlg.SetSize(Size)
                dlg.CenterOnParent()

                # publication info
                step = 1
                dlg.Update(step,"Exporting overall section")
                WriteCIFitem('\ndata_'+self.CIFname+'_publ')
                WriteAudit()
                WriteCIFitem('_pd_block_id',
                             str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                             str(self.shortauthorname) + "|Overall")
                writeCIFtemplate(self.OverallParms['Controls'],'publ') #insert the publication template
                # ``template_publ.cif`` or a modified version
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
                    j = histblk['hId']
                    datablockidDict[hist] = (str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                                             str(self.shortauthorname) + "|" +
                                             instnam + "_hist_"+str(j))
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
                    step += 1
                    dlg.Update(step,"Exporting phase "+phasenam+' (#'+str(j+1)+')')
                    i = self.Phases[phasenam]['pId']
                    WriteCIFitem('\ndata_'+self.CIFname+"_phase_"+str(i))
                    WriteCIFitem('# Information for phase '+str(i))
                    WriteCIFitem('_pd_block_id',datablockidDict[phasenam])
                    # report the phase
                    writeCIFtemplate(self.Phases[phasenam]['General'],'phase',phasenam) # write phase template
                    WritePhaseInfo(phasenam)
                    # preferred orientation
                    if self.ifPWDR:
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
                        step += 1
                        dlg.Update(step,"Exporting "+hist.strip())
                        WriteCIFitem('\ndata_'+self.CIFname+"_pwd_"+str(i))
                        #instnam = histblk["Sample Parameters"]['InstrName']
                        # report instrumental profile terms
                        WriteCIFitem('_pd_proc_ls_profile_function',
                            FormatInstProfile(histblk["Instrument Parameters"],histblk['hId']))
                        WriteCIFitem('# Information for histogram '+str(i)+': '+hist)
                        WriteCIFitem('_pd_block_id',datablockidDict[hist])
                        histprm = self.Histograms[hist]["Sample Parameters"]
                        writeCIFtemplate(histprm,'powder',histprm['InstrName']) # powder template
                        WritePowderData(hist)
                for i in sorted(self.xtalDict.keys()):
                    hist = self.xtalDict[i]
                    histblk = self.Histograms[hist]
                    if hist.startswith("HKLF"): 
                        step += 1
                        dlg.Update(step,"Exporting "+hist.strip())
                        WriteCIFitem('\ndata_'+self.CIFname+"_sx_"+str(i))
                        #instnam = histblk["Instrument Parameters"][0]['InstrName']
                        WriteCIFitem('# Information for histogram '+str(i)+': '+hist)
                        WriteCIFitem('_pd_block_id',datablockidDict[hist])
                        histprm = self.Histograms[hist]["Instrument Parameters"][0]
                        writeCIFtemplate(histprm,'single',histprm['InstrName']) # single crystal template
                        WriteSingleXtalData(hist)

            except Exception:
                import traceback
                print(traceback.format_exc())
                self.G2frame.ErrorDialog('Exception',
                                         'Error occurred in CIF creation. '+
                                         'See full error message in console output ')
            finally:
                dlg.Destroy()

        WriteCIFitem('#--' + 15*'eof--' + '#')
        self.CloseFile()
        print("...export completed")
        print('file '+str(self.fullpath))
        # end of CIF export

class ExportPhaseCIF(ExportCIF):
    '''Used to create a simple CIF of at most one phase. Uses exact same code as
    :class:`ExportCIF` except that `self.mode` is set to "simple" in `self.InitExport`.
    Shows up in menu as Quick CIF.

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        G2IO.ExportBaseclass.__init__(self,
            G2frame=G2frame,
            formatName = 'Quick CIF',
            extension='.cif',
            longFormatName = 'Export one phase in CIF'
            )
        self.exporttype = ['phase']
        # CIF-specific items
        self.author = ''
        self.mode = 'simple'

class ExportDataCIF(ExportCIF):
    '''Used to create a simple CIF containing diffraction data only. Uses exact same code as
    :class:`ExportCIF` except that `self.mode` is set to "simple" and `self.currentExportType`
    is set to "single" or "powder" in `self.InitExport`. Shows up in menus as Data-only CIF.

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        G2IO.ExportBaseclass.__init__(self,
            G2frame=G2frame,
            formatName = 'Data-only CIF',
            extension='.cif',
            longFormatName = 'Export data as CIF'
            )
        self.exporttype = ['single','powder']
        # CIF-specific items
        self.author = ''
        self.mode = 'simple'

#===============================================================================
# misc CIF utilities
#===============================================================================
def PickleCIFdict(fil):
    '''Loads a CIF dictionary, cherry picks out the items needed
    by local code and sticks them into a python dict and writes
    that dict out as a cPickle file for later reuse.
    If the write fails a warning message is printed,
    but no exception occurs.

    :param str fil: file name of CIF dictionary, will usually end
      in .dic
    :returns: the dict with the definitions  
    '''
    import CifFile as cif # PyCifRW from James Hester
    cifdic = {}
    try:
        fp = open(fil,'r')             # patch: open file to avoid windows bug
        dictobj = cif.CifDic(fp)
        fp.close()
    except IOError:
        dictobj = cif.CifDic(fil)
    if DEBUG: print('loaded '+str(fil))
    for item in dictobj.keys():
        cifdic[item] = {}
        for j in (
            '_definition','_type',
            '_enumeration',
            '_enumeration_detail',
            '_enumeration_range'):
            if dictobj[item].get(j):
                cifdic[item][j] = dictobj[item][j]
    try:
        fil = os.path.splitext(fil)[0]+'.cpickle'
        fp = open(fil,'w')
        cPickle.dump(cifdic,fp)
        fp.close()
        if DEBUG: print('wrote '+str(fil))
    except:
        print ('Unable to write '+str(fil))
    return cifdic

def LoadCIFdic():
    '''Create a composite core+powder CIF lookup dict containing
    information about all items in the CIF dictionaries, loading
    pickled files if possible. The routine looks for files
    named cif_core.cpickle and cif_pd.cpickle in every
    directory in the path and if they are not found, files
    cif_core.dic and/or cif_pd.dic are read.

    :returns: the dict with the definitions  
    '''
    cifdic = {}
    for ftyp in "cif_core","cif_pd":
        for loc in sys.path:
            fil = os.path.join(loc,ftyp+".cpickle")
            if not os.path.exists(fil): continue
            fp = open(fil,'r')
            try:
                cifdic.update(cPickle.load(fp))
                if DEBUG: print('reloaded '+str(fil))
                break
            finally:
                fp.close()
        else:
            for loc in sys.path:
                fil = os.path.join(loc,ftyp+".dic")
                if not os.path.exists(fil): continue
                #try:
                if True:
                    cifdic.update(PickleCIFdict(fil))
                    break
                #except:
                #    pass
            else:
                print('Could not load '+ftyp+' dictionary')
    return cifdic

class CIFdefHelp(wx.Button):
    '''Create a help button that displays help information on
    the current data item

    :param parent: the panel which will be the parent of the button
    :param str msg: the help text to be displayed
    :param wx.Dialog helpwin: Frame for CIF editing dialog
    :param wx.TextCtrl helptxt: TextCtrl where help text is placed
    '''
    def __init__(self,parent,msg,helpwin,helptxt):
        wx.Button.__init__(self,parent,wx.ID_HELP)
        self.Bind(wx.EVT_BUTTON,self._onPress)
        self.msg=msg
        self.parent = parent
        #self.helpwin = self.parent.helpwin
        self.helpwin = helpwin
        self.helptxt = helptxt
    def _onPress(self,event):
        'Respond to a button press by displaying the requested text'
        try:
            #helptxt = self.helptxt
            ow,oh = self.helptxt.GetSize()
            self.helptxt.SetLabel(self.msg)
            w,h = self.helptxt.GetSize()
            if h > oh:
                self.helpwin.GetSizer().Fit(self.helpwin)
        except: # error posting help, ignore
            return

def CIF2dict(cf):
    '''copy the contents of a CIF out from a PyCifRW block object
    into a dict

    :returns: cifblk, loopstructure where cifblk is a dict with
      CIF items and loopstructure is a list of lists that defines
      which items are in which loops. 
    '''
    blk = cf.keys()[0] # assume templates are a single CIF block, use the 1st
    loopstructure = cf[blk].loopnames()[:] # copy over the list of loop contents
    dblk = {}
    for item in cf[blk].keys(): # make a copy of all the items in the block
        dblk[item] = cf[blk][item]
    return dblk,loopstructure

def dict2CIF(dblk,loopstructure,blockname='Template'):
    '''Create a PyCifRW CIF object containing a single CIF
    block object from a dict and loop structure list.

    :param dblk: a dict containing values for each CIF item
    :param list loopstructure: a list of lists containing the contents of
      each loop, as an example::

         [ ["_a","_b"], ["_c"], ["_d_1","_d_2","_d_3"]]

      this describes a CIF with this type of structure::

        loop_ _a _b <a1> <b1> <a2> ...
        loop_ _c <c1> <c2>...
        loop _d_1 _d_2 _d_3 ...

      Note that the values for each looped CIF item, such as _a,
      are contained in a list, for example as cifblk["_a"]

    :param str blockname: an optional name for the CIF block.
      Defaults to 'Template'

    :returns: the newly created PyCifRW CIF object 
    '''

    import CifFile as cif # PyCifRW from James Hester
    # compile a 'list' of items in loops
    loopnames = set()
    for i in loopstructure:
        loopnames |= set(i)
    # create a new block
    newblk = cif.CifBlock()
    # add the looped items
    for keys in loopstructure:
        vals = []
        for key in keys:
            vals.append(dblk[key])
        newblk.AddCifItem(([keys],[vals]))
    # add the non-looped items
    for item in dblk:
        if item in loopnames: continue
        newblk[item] = dblk[item]
    # create a CIF and add the block
    newcf = cif.CifFile()
    newcf[blockname] = newblk    
    return newcf


class EditCIFtemplate(wx.Dialog):
    '''Create a dialog for editing a CIF template. The edited information is
    placed in cifblk. If the CIF is saved as a file, the name of that file
    is saved as ``self.newfile``.
    
    :param wx.Frame parent: parent frame or None
    :param cifblk: dict or PyCifRW block containing values for each CIF item
    :param list loopstructure: a list of lists containing the contents of
      each loop, as an example::

         [ ["_a","_b"], ["_c"], ["_d_1","_d_2","_d_3"]]

      this describes a CIF with this type of structure::

        loop_ _a _b <a1> <b1> <a2> ...
        loop_ _c <c1> <c2>...
        loop _d_1 _d_2 _d_3 ...

      Note that the values for each looped CIF item, such as _a,
      are contained in a list, for example as cifblk["_a"]
      
    :param str defaultname: specifies the default file name to be used for 
      saving the CIF.
    '''
    def __init__(self,parent,cifblk,loopstructure,defaultname):
        OKbuttons = []
        self.cifblk = cifblk
        self.loopstructure = loopstructure
        self.newfile = None
        self.defaultname = defaultname        
        global CIFdic  # once this is loaded, keep it around
        if CIFdic is None:
            CIFdic = LoadCIFdic()
        wx.Dialog.__init__(self,parent,style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        # define widgets that will be needed during panel creation
        self.helptxt = wx.StaticText(self,wx.ID_ANY,"")
        savebtn = wx.Button(self, wx.ID_CLOSE, "Save as template")
        OKbuttons.append(savebtn)
        savebtn.Bind(wx.EVT_BUTTON,self._onSave)
        OKbtn = wx.Button(self, wx.ID_OK, "Use")
        OKbtn.SetDefault()
        OKbuttons.append(OKbtn)

        self.SetTitle('Edit items in CIF template')
        vbox = wx.BoxSizer(wx.VERTICAL)
        cpnl = EditCIFpanel(self,cifblk,loopstructure,CIFdic,OKbuttons,size=(300,300))
        vbox.Add(cpnl, 1, wx.ALIGN_LEFT|wx.ALL|wx.EXPAND, 0)
        G2gd.HorizontalLine(vbox,self)
        vbox.Add(self.helptxt, 0, wx.EXPAND|wx.ALL, 5)
        G2gd.HorizontalLine(vbox,self)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.Add(btn,0,wx.ALIGN_CENTER|wx.ALL)
        btnsizer.Add(savebtn,0,wx.ALIGN_CENTER|wx.ALL)
        btnsizer.Add(OKbtn,0,wx.ALIGN_CENTER|wx.ALL)
        vbox.Add(btnsizer, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        self.SetSizer(vbox)
        vbox.Fit(self)
    def Post(self):
        '''Display the dialog
        
        :returns: True unless Cancel has been pressed.
        '''
        return (self.ShowModal() == wx.ID_OK)
    def _onSave(self,event):
        'Save CIF entries in a template file'
        dlg = wx.FileDialog(
            self, message="Save as CIF template",
            defaultDir=os.getcwd(),
            defaultFile=self.defaultname,
            wildcard="CIF (*.cif)|*.cif",
            style=wx.SAVE | wx.CHANGE_DIR
            )
        val = (dlg.ShowModal() == wx.ID_OK)
        fil = dlg.GetPath()
        dlg.Destroy()
        if val: # ignore a Cancel button
            fil = os.path.splitext(fil)[0]+'.cif' # force extension
            fp = open(fil,'w')
            newcf = dict2CIF(self.cifblk,self.loopstructure)
            fp.write(newcf.WriteOut())
            fp.close()
            self.newfile = fil
            self.EndModal(wx.ID_OK)

class EditCIFpanel(wxscroll.ScrolledPanel):
    '''Creates a scrolled panel for editing CIF template items

    :param wx.Frame parent: parent frame where panel will be placed
    :param cifblk: dict or PyCifRW block containing values for each CIF item
    :param list loopstructure: a list of lists containing the contents of
      each loop, as an example::

         [ ["_a","_b"], ["_c"], ["_d_1","_d_2","_d_3"]]

      this describes a CIF with this type of structure::

        loop_ _a _b <a1> <b1> <a2> ...
        loop_ _c <c1> <c2>...
        loop _d_1 _d_2 _d_3 ...

      Note that the values for each looped CIF item, such as _a,
      are contained in a list, for example as cifblk["_a"]

    :param dict cifdic: optional CIF dictionary definitions
    :param list OKbuttons: A list of wx.Button objects that should
      be disabled when information in the CIF is invalid
    :param (other): optional keyword parameters for wx.ScrolledPanel
    '''
    def __init__(self, parent, cifblk, loopstructure, cifdic={}, OKbuttons=[], **kw):
        self.parent = parent
        wxscroll.ScrolledPanel.__init__(self, parent, wx.ID_ANY, **kw)
        self.vbox = None
        self.AddDict = None
        self.cifdic = cifdic
        self.cifblk = cifblk
        self.loops = loopstructure
        self.parent = parent
        self.LayoutCalled = False
        self.parentOKbuttons = OKbuttons
        self.ValidatedControlsList = []
        self._fill()
    def _fill(self):
        'Fill the scrolled panel with widgets for each CIF item'
        wx.BeginBusyCursor()
        self.AddDict = {}
        self.ValidatedControlsList = []
        # delete any only contents
        if self.vbox:
            self.vbox.DeleteWindows()
            self.vbox = None
            self.Update()
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox = vbox
        # compile a 'list' of items in loops
        loopnames = set()
        for i in self.loops:
            loopnames |= set(i)
        # post the looped CIF items
        for lnum,lp in enumerate(self.loops):
            hbox = wx.BoxSizer(wx.HORIZONTAL)
            hbox.Add(wx.StaticText(self,wx.ID_ANY,'Loop '+str(lnum+1)))
            vbox.Add(hbox)
            but = wx.Button(self,wx.ID_ANY,"Add row")
            self.AddDict[but]=lnum
            
            hbox.Add(but)            
            but.Bind(wx.EVT_BUTTON,self.OnAddRow)
            fbox = wx.GridBagSizer(0, 0)
            vbox.Add(fbox)
            rows = 0
            for i,item in enumerate(lp):
                txt = wx.StaticText(self,wx.ID_ANY,item+"  ")
                fbox.Add(txt,(0,i+1))
                # if self.cifdic.get(item):
                #     df = self.cifdic[item].get('_definition')
                #     if df:
                #         txt.SetToolTipString(G2IO.trim(df))
                #         but = CIFdefHelp(self,
                #                          "Definition for "+item+":\n\n"+G2IO.trim(df),
                #                          self.parent,
                #                          self.parent.helptxt)
                #         fbox.Add(but,(1,i+1),flag=wx.ALIGN_CENTER)
                for j,val in enumerate(self.cifblk[item]):
                    ent = self.CIFEntryWidget(self.cifblk[item],j,item)
                    #fbox.Add(ent,(j+2,i+1),flag=wx.EXPAND|wx.ALL)
                    fbox.Add(ent,(j+1,i+1),flag=wx.EXPAND|wx.ALL)
                if self.cifdic.get(item):
                    df = self.cifdic[item].get('_definition')
                    if df:
                        txt.SetToolTipString(G2IO.trim(df))
                        but = CIFdefHelp(self,
                                         "Definition for "+item+":\n\n"+G2IO.trim(df),
                                         self.parent,
                                         self.parent.helptxt)
                        fbox.Add(but,(j+2,i+1),flag=wx.ALIGN_CENTER)
                rows = max(rows,len(self.cifblk[item]))
            for i in range(rows):
                txt = wx.StaticText(self,wx.ID_ANY,str(i+1))
                fbox.Add(txt,(i+2,0))
            line = wx.StaticLine(self,wx.ID_ANY, size=(-1,3), style=wx.LI_HORIZONTAL)
            vbox.Add(line, 0, wx.EXPAND|wx.ALIGN_CENTER|wx.ALL, 10)
                
        # post the non-looped CIF items
        for item in sorted(self.cifblk.keys()):
            if item not in loopnames:
                hbox = wx.BoxSizer(wx.HORIZONTAL)
                vbox.Add(hbox)
                txt = wx.StaticText(self,wx.ID_ANY,item)
                hbox.Add(txt)
                ent = self.CIFEntryWidget(self.cifblk,item,item)
                hbox.Add(ent)
                if self.cifdic.get(item):
                    df = self.cifdic[item].get('_definition')
                    if df:
                        txt.SetToolTipString(G2IO.trim(df))
                        but = CIFdefHelp(self,
                                         "Definition for "+item+":\n\n"+G2IO.trim(df),
                                         self.parent,
                                         self.parent.helptxt)
                        hbox.Add(but,0,wx.ALL,2)
        self.SetSizer(vbox)
        #vbox.Fit(self.parent)
        self.SetAutoLayout(1)
        self.SetupScrolling()
        self.Bind(rw.EVT_RW_LAYOUT_NEEDED, self.OnLayoutNeeded)
        self.Layout()
        wx.EndBusyCursor()
    def OnLayoutNeeded(self,event):
        '''Called when an update of the panel layout is needed. Calls
        self.DoLayout after the current operations are complete using
        CallAfter. This is called only once, according to flag
        self.LayoutCalled, which is cleared in self.DoLayout. 
        '''
        if self.LayoutCalled: return # call already queued
        wx.CallAfter(self.DoLayout) # queue a call
        self.LayoutCalled = True
    def DoLayout(self):
        '''Update the Layout and scroll bars for the Panel. Clears
        self.LayoutCalled so that next change to panel can
        request a new update
        '''
        wx.BeginBusyCursor()
        self.Layout()
        self.SetupScrolling()
        wx.EndBusyCursor()
        self.LayoutCalled = False
    def OnAddRow(self,event):
        'add a row to a loop'
        lnum = self.AddDict.get(event.GetEventObject())
        if lnum is None: return
        for item in self.loops[lnum]:
            self.cifblk[item].append('?')
        self._fill()

    def ControlOKButton(self,setvalue):
        '''Enable or Disable the OK button(s) for the dialog. Note that this is
        passed into the ValidatedTxtCtrl for use by validators.

        :param bool setvalue: if True, all entries in the dialog are
          checked for validity. The first invalid control triggers
          disabling of buttons.
          If False then the OK button(s) are disabled with no checking
          of the invalid flag for each control.
        '''
        if setvalue: # turn button on, do only if all controls show as valid
            for ctrl in self.ValidatedControlsList:
                if ctrl.invalid:
                    for btn in self.parentOKbuttons:
                        btn.Disable()
                    return
            else:
                for btn in self.parentOKbuttons:
                    btn.Enable()
        else:
            for btn in self.parentOKbuttons:
                btn.Disable()
        
    def CIFEntryWidget(self,dct,item,dataname):
        '''Create an entry widget for a CIF item. Use a validated entry for numb values
        where int is required when limits are integers and floats otherwise.
        At present this does not allow entry of the special CIF values of "." and "?" for
        numerical values and highlights them as invalid.
        Use a selection widget when there are specific enumerated values for a string.        
        '''
        if self.cifdic.get(dataname):
            if self.cifdic[dataname].get('_enumeration'):
                values = ['?']+self.cifdic[dataname]['_enumeration']
                choices = ['undefined']
                for i in self.cifdic[dataname].get('_enumeration_detail',values):
                    choices.append(G2IO.trim(i))
                ent = G2gd.EnumSelector(self, dct, item, choices, values, size=(200, -1))
                return ent
            if self.cifdic[dataname].get('_type') == 'numb':
                mn = None
                mx = None
                hint = int
                if self.cifdic[dataname].get('_enumeration_range'):
                    rng = self.cifdic[dataname]['_enumeration_range'].split(':')
                    if '.' in rng[0] or '.' in rng[1]: hint = float
                    if rng[0]: mn = hint(rng[0])
                    if rng[1]: mx = hint(rng[1])
                    ent = G2ctrls.ValidatedTxtCtrl(
                        self,dct,item,typeHint=hint,min=mn,max=mx,
                        CIFinput=True,
                        OKcontrol=self.ControlOKButton)
                    self.ValidatedControlsList.append(ent)
                    return ent
        rw1 = rw.ResizeWidget(self)
        ent = G2ctrls.ValidatedTxtCtrl(
            rw1,dct,item,size=(100, 20),
            style=wx.TE_MULTILINE|wx.TE_PROCESS_ENTER,
            CIFinput=True,
            OKcontrol=self.ControlOKButton)
        self.ValidatedControlsList.append(ent)
        return rw1

class CIFtemplateSelect(wx.BoxSizer):
    '''Create a set of buttons to show, select and edit a CIF template
    
    :param frame: wx.Frame object of parent
    :param panel: wx.Panel object where widgets should be placed
    :param str tmplate: one of 'publ', 'phase', or 'instrument' to determine
      the type of template
    :param dict G2dict: GSAS-II dict where CIF should be placed. The key
      "CIF_template" will be used to store either a list or a string.
      If a list, it will contain a dict and a list defining loops. If
      an str, it will contain a file name.    
    :param function repaint: reference to a routine to be called to repaint
      the frame after a change has been made
    :param str title: A line of text to show at the top of the window
    :param str defaultname: specifies the default file name to be used for 
      saving the CIF.
    '''
    def __init__(self,frame,panel,tmplate,G2dict, repaint, title, defaultname=''):
        wx.BoxSizer.__init__(self,wx.VERTICAL)
        self.cifdefs = frame
        self.dict = G2dict
        self.repaint = repaint
        templateDefName = 'template_'+tmplate+'.cif'
        self.CIF = G2dict.get("CIF_template")
        if defaultname:
            self.defaultname = defaultname.encode('ascii','replace').strip().replace(' ','_')
            self.defaultname = re.sub(r'[^a-zA-Z0-9_-]','',self.defaultname)
            self.defaultname = tmplate + "_" + self.defaultname + ".cif"
        else:
            self.defaultname = ''
            
        txt = wx.StaticText(panel,wx.ID_ANY,title)
        self.Add(txt,0,wx.ALIGN_CENTER)
        # change font on title
        txtfnt = txt.GetFont()
        txtfnt.SetWeight(wx.BOLD)
        txtfnt.SetPointSize(2+txtfnt.GetPointSize())
        txt.SetFont(txtfnt)
        self.Add((-1,3))

        if not self.CIF: # empty or None
            for pth in [os.getcwd()]+sys.path:
                fil = os.path.join(pth,self.defaultname)
                if os.path.exists(fil) and self.defaultname:
                    self.CIF = fil
                    CIFtxt = "Template: "+self.defaultname
                    break
            else:
                for pth in sys.path:
                    fil = os.path.join(pth,templateDefName)
                    if os.path.exists(fil):
                        self.CIF = fil
                        CIFtxt = "Template: "+templateDefName
                        break
                else:
                    print("Default CIF template "+self.defaultname+' not found in path!')
                    self.CIF = None
                    CIFtxt = "none! (No template found)"
        elif type(self.CIF) is not list and type(self.CIF) is not tuple:
            if not os.path.exists(self.CIF):
                print("Error: template file has disappeared: "+self.CIF)
                self.CIF = None
                CIFtxt = "none! (file not found)"
            else:
                if len(self.CIF) < 50:
                    CIFtxt = "File: "+self.CIF
                else:
                    CIFtxt = "File: ..."+self.CIF[-50:]
        else:
            CIFtxt = "Template is customized"
        # show template source
        self.Add(wx.StaticText(panel,wx.ID_ANY,CIFtxt))
        # show str, button to select file; button to edit (if CIF defined)
        but = wx.Button(panel,wx.ID_ANY,"Select Template File")
        but.Bind(wx.EVT_BUTTON,self._onGetTemplateFile)
        hbox =  wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(but,0,0,2)
        but = wx.Button(panel,wx.ID_ANY,"Edit Template")
        but.Bind(wx.EVT_BUTTON,self._onEditTemplateContents)
        if self.CIF is None: but.Disable() # nothing to edit!
        hbox.Add(but,0,0,2)
        self.Add(hbox)
    def _onGetTemplateFile(self,event):
        'select a template file'
        dlg = wx.FileDialog(
            self.cifdefs, message="Read CIF template file",
            defaultDir=os.getcwd(),
            defaultFile=self.defaultname,
            wildcard="CIF (*.cif)|*.cif",
            style=wx.OPEN | wx.CHANGE_DIR
            )
        ret = dlg.ShowModal()
        fil = dlg.GetPath()
        dlg.Destroy()
        if ret == wx.ID_OK:
            try:
                cf = G2IO.ReadCIF(fil)
                if len(cf.keys()) == 0:
                    raise Exception,"No CIF data_ blocks found"
                if len(cf.keys()) != 1:
                    raise Exception, 'Error, CIF Template has more than one block: '+fil
                self.dict["CIF_template"] = fil
            except Exception as err:
                print('\nError reading CIF: '+fil)
                dlg = wx.MessageDialog(self.cifdefs,
                                   'Error reading CIF '+fil,
                                   'Error in CIF file',
                                   wx.OK)
                dlg.ShowModal()
                dlg.Destroy()
                print(err.message)
                return
            self.repaint() #EditCIFDefaults()

    def _onEditTemplateContents(self,event):
        'Called to edit the contents of a CIF template'
        if type(self.CIF) is list or  type(self.CIF) is tuple:
            dblk,loopstructure = copy.deepcopy(self.CIF) # don't modify original
        else:
            cf = G2IO.ReadCIF(self.CIF)
            dblk,loopstructure = CIF2dict(cf)
        dlg = EditCIFtemplate(self.cifdefs,dblk,loopstructure,self.defaultname)
        val = dlg.Post()
        if val:
            if dlg.newfile: # results saved in file
                self.dict["CIF_template"] = dlg.newfile
            else:
                self.dict["CIF_template"] = [dlg.cifblk,dlg.loopstructure]
            self.repaint() #EditCIFDefaults() # note that this does a dlg.Destroy()
        else:
            dlg.Destroy()        

#===============================================================================
# end of misc CIF utilities
#===============================================================================
