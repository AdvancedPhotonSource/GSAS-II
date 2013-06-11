'''Development code to export a GSAS-II project as a CIF
The heavy lifting is done in method export
'''

import datetime as dt
import os.path
import GSASIIIO as G2IO
import GSASIIgrid as G2gd
#reload(G2gd)
class ExportCIF(G2IO.ExportBaseclass):
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'full CIF',
            longFormatName = 'Export project as complete CIF'
            )
        self.author = ''

    def export(self,mode='full'):
        '''Export a CIF

        :param str mode: "full" (default) to create a complete CIF of project,
          "simple" for a simple CIF with only coordinates
        '''
        def WriteCIFitem(name,value=''):
            if value:
                if "\n" in value or len(value)> 70:
                    if name.strip(): print name
                    print '; '+value
                    print '; '
                elif " " in value:
                    if len(name)+len(value) > 65:
                        print name,'\n   ','"' + str(value) + '"'
                    else:
                        print name,'  ','"' + str(value) + '"'
                else:
                    if len(name)+len(value) > 65:
                        print name,'\n   ',value
                    else:
                        print name,'  ',value
            else:
                print name

        def WriteAudit():
            WriteCIFitem('_audit_creation_method',
                         'created in GSAS-II')
            WriteCIFitem('_audit_creation_date',self.CIFdate)
            if self.author:
                WriteCIFitem('_audit_author_name',self.author)
            WriteCIFitem('_audit_update_record',
                         self.CIFdate+'  Initial software-generated CIF')

        def WriteOverall():
            '''TODO: Write out overall refinement information
            '''
            WriteCIFitem('_pd_proc_info_datetime', self.CIFdate)
            WriteCIFitem('_pd_calc_method', 'Rietveld Refinement')
            #WriteCIFitem('_refine_ls_shift/su_max',DAT1)
            #WriteCIFitem('_refine_ls_shift/su_mean',DAT2)
            WriteCIFitem('_computing_structure_refinement','GSAS-II')
            #WriteCIFitem('_refine_ls_number_parameters',DAT1)
            #WriteCIFitem('_refine_ls_goodness_of_fit_all',DAT2)
            #WriteCIFitem('_refine_ls_number_restraints',TEXT(1:7))
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
            if self.npowder > 1:
                WriteCIFitem('\n# OVERALL POWDER R-FACTORS')
                #WriteCIFitem('_pd_proc_ls_prof_R_factor',TEXT(11:20))
                #WriteCIFitem('_pd_proc_ls_prof_wR_factor',TEXT(1:10))
            WriteCIFitem('_refine_ls_matrix_type','full')
            #WriteCIFitem('_refine_ls_matrix_type','userblocks')

        def WritePubTemplate():
            '''TODO: insert the publication template ``template_publ.cif`` or some modified
            version for this project. Store this in the GPX file?
            '''
            print "TODO: publication info goes here"

        def WritePhaseTemplate():
            '''TODO: insert the phase template ``template_phase.cif`` or some modified
            version for this project
            '''
            print "TODO: phase template info goes here"

        def WritePowderTemplate():
            '''TODO: insert the phase template ``template_instrument.cif`` or some modified
            version for this project
            '''
            print "TODO: powder histogram template info goes here"

        def WriteSnglXtalTemplate():
            '''TODO: insert the single-crystal histogram template 
            for this project
            '''
            print "TODO: single-crystal histogram template info goes here"

        def WritePhaseInfo(phasenam):
            # see writepha.for
            print 'TODO: phase info for',phasenam,'goes here'
            # THINK: how to select publication flags for distances/angles?
            phasedict = self.GroupedParms['Phases'][phasenam] # pointer to current phase info            
            WriteCIFitem('_pd_phase_name', phasenam)
            text = '?'
            for lbl in 'a','b','c':
                WriteCIFitem('_cell_length_'+lbl,text)
            for lbl in 'alpha','beta ','gamma':
                WriteCIFitem('_cell_angle_'+lbl,text)

            WriteCIFitem('_cell_volume',text)

            #print phasedict['Histograms']
#            for key in phasedict['General']:
#                print key,'=',phasedict['General'][key]
#            print phasedict['pId']

            #NSYS = (1,2,3,4,4,5,5,5,5,5,6,6,7,7)
            #SYST = ('','triclinic','monoclinic','orthorhombic','tetragonal','trigonal ','hexagonal','cubic')

            WriteCIFitem('_symmetry_cell_setting',
                         phasedict['General']['SGData']['SGSys'])

            spacegroup = phasedict['General']['SGData']['SpGrp'].strip()
            # regularize capitalization and remove trailing H/R
            spacegroup = spacegroup[0].upper() + spacegroup[1:].lower().rstrip('rh ')
            WriteCIFitem('_symmetry_space_group_name_H-M',spacegroup)
            # generate symmetry operations including centering and center of symmetry
            WriteCIFitem('loop_ _symmetry_equiv_pos_site_id _symmetry_equiv_pos_as_xyz')

            # loop over histogram(s) used in this phase
            if not oneblock:
                # pointers to histograms used in this phase
                histlist = []
                for hist in self.GroupedParms['Phases'][phasenam]['Histograms']:
                    if self.GroupedParms['Phases'][phasenam]['Histograms'][hist]['Use']:
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
                    for hist in histlist:
                        WriteCIFitem('',hist)
                        # TODO: sample/histogram profile information and other "HAP" info should be
                        # reported in this loop, such as
                        # _pd_proc_ls_pref_orient_corr 
                        # _pd_proc_ls_profile_function
                        #txt = ' Spherical Harmonic ODF' + '\n spherical harmonic order='+'?'
                        #IF ( ISAMSYM.EQ.1 ) THEN
                        #txt += '\n No sample symmetry'
                #ELSE IF ( ISAMSYM.EQ.2 ) THEN
                #txt += ' The sample symmetry is: 2/m (shear texture)'
                #ELSE IF ( ISAMSYM.EQ.3 ) THEN
                #txt += ' The sample symmetry is: mmm (rolling texture)'
                #ELSE IF ( ISAMSYM.EQ.0 ) THEN
                #txt += ' The sample symmetry is: cylindrical (fiber texture)'
                #WriteCIFitem('_pd_proc_ls_pref_orient_corr',txt)
            else:
                # report all profile information here (include instrumental)
                # _pd_proc_ls_pref_orient_corr 
                # _pd_proc_ls_profile_function 
                pass
            
            WriteCIFitem('\n# ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS')
      
            WriteCIFitem('loop_ '+
                         '\n\t_atom_site_type_symbol'+
                         '\n\t_atom_site_label'+
                         '\n\t_atom_site_fract_x'+
                         '\n\t_atom_site_fract_y'+
                         '\n\t_atom_site_fract_z'+
                         '\n\t_atom_site_occupancy'+
                         '\n\t_atom_site_thermal_displace_type'+
                         '\n\t_atom_site_U_iso_or_equiv'+
                         '\n\t_atom_site_symmetry_multiplicity')

            WriteCIFitem('loop_' + '\n\t_atom_site_aniso_label' + 
                         '\n\t_atom_site_aniso_U_11' + '\n\t_atom_site_aniso_U_12' +
                         '\n\t_atom_site_aniso_U_13' + '\n\t_atom_site_aniso_U_22' +
                         '\n\t_atom_site_aniso_U_23' + '\n\t_atom_site_aniso_U_33')

            # process the chemical formula: pick a Z value & generate molecular weight
            # find the maximum possible Z value

            # order the elements in "Hill" order: C, H, D, T & alphabetical or alphabetical
            if not oneblock: # in single block, this is combined with the scattering factors
                WriteCIFitem('loop_  _atom_type_symbol _atom_type_number_in_cell')

            WriteCIFitem('# If you change Z, be sure to change all 3 of the following')
            WriteCIFitem( '_chemical_formula_sum',text)
            #WRITE(TEXT,'(F15.2)') ATMASS
            WriteCIFitem( '_chemical_formula_weight',text)
            #WRITE(TEXT,'(I4)') Z
            WriteCIFitem( '_cell_formula_units_Z',text)

            #C now loop over interatomic distances for this phase
            WriteCIFitem('\n# MOLECULAR GEOMETRY')
            WriteCIFitem('loop_' + 
                         '\n\t_geom_bond_atom_site_label_1' +
                         '\n\t_geom_bond_atom_site_label_2' + 
                         '\n\t_geom_bond_distance' + 
                         '\n\t_geom_bond_site_symmetry_1' + 
                         '\n\t_geom_bond_site_symmetry_2' + 
                         '\n\t_geom_bond_publ_flag')

            #C now loop over interatomic angles for this phase
            WriteCIFitem('loop_' + 
                         '\n\t_geom_angle_atom_site_label_1' + 
                         '\n\t_geom_angle_atom_site_label_2' + 
                         '\n\t_geom_angle_atom_site_label_3' + 
                         '\n\t_geom_angle' + 
                         '\n\t_geom_angle_site_symmetry_1' +
                         '\n\t_geom_angle_site_symmetry_2' + 
                         '\n\t_geom_angle_site_symmetry_3' + 
                         '\n\t_geom_angle_publ_flag')

        def WritePowderData(histlbl):
            text = '?'
            histblk = self.GroupedParms['PWDR'][histlbl]
            print 'TODO: powder here data for',histblk["Sample Parameters"]['InstrName']
            # see wrpowdhist.for & wrreflist.for
            
            refprx = '_refln.' # mm
            refprx = '_refln_' # normal

            if not oneblock:
                if not phasebyhistDict.get(histlbl):
                    WriteCIFitem('\n# No phases associated with this data set')
                else:
                    WriteCIFitem('\n# PHASE TABLE')
                    WriteCIFitem('loop_' +
                                 '\n\t_pd_phase_id' + 
                                 '\n\t_pd_phase_block_id' + 
                                 '\n\t_pd_phase_mass_%')
                    for phasenam in phasebyhistDict.get(histlbl):
                        pass

            WriteCIFitem('\n# SCATTERING FACTOR INFO')
            WriteCIFitem('_diffrn_radiation_wavelength' ,text)
            #WriteCIFitem('_diffrn_radiation_type',text)
            #C always assume Ka1 & Ka2 if two wavelengths are present
            #WriteCIFitem('loop_' + 
            #             '\n\t_diffrn_radiation_wavelength' +
            #             '\n\t_diffrn_radiation_wavelength_wt' + 
            #             '\n\t_diffrn_radiation_type' + 
            #             '\n\t_diffrn_radiation_wavelength_id')
            #WRITE LAM1,1.0,'K\\a~1~',1, LAM2,ratio,'K\\a~2~',2

            WriteCIFitem('_pd_proc_ls_prof_R_factor','?')
            WriteCIFitem('_pd_proc_ls_prof_wR_factor','?')
            WriteCIFitem('_pd_proc_ls_prof_wR_expected','?')
            WriteCIFitem('_refine_ls_R_Fsqd_factor','?')
            
            #WriteCIFitem('_pd_meas_2theta_fixed',text)
            WriteCIFitem('_diffrn_radiation_probe','x-ray')
            WriteCIFitem('_diffrn_radiation_probe','neutron')
            WriteCIFitem('_diffrn_radiation_polarisn_ratio','?')
            
            WriteCIFitem('loop_  _atom_type_symbol')
            if oneblock:
                WriteCIFitem('       _atom_type_number_in_cell')
            #IF (HTYP(2:2) .eq. 'X' .AND. HTYP(3:3) .ne. 'E') THEN
            WriteCIFitem('      _atom_type_scat_dispersion_real')
            WriteCIFitem('      _atom_type_scat_dispersion_imag')
            for lbl in ('a1','a2','a3', 'a4', 'b1', 'b2', 'b3', 'b4', 'c'):
                WriteCIFitem('      _atom_type_scat_Cromer_Mann_'+lbl)
            #ELSEIF (HTYP(2:2) .eq. 'N') THEN
            WriteCIFitem('      _atom_type_scat_length_neutron')
            #ENDIF
            WriteCIFitem('      _atom_type_scat_source')

            #C document the background function used
            WriteCIFitem('_pd_proc_ls_background_function','?')

            WriteCIFitem('_exptl_absorpt_process_details','?')
            WriteCIFitem('_exptl_absorpt_correction_T_min','?')
            WriteCIFitem('_exptl_absorpt_correction_T_max','?')
            #C extinction
            #WRITE(IUCIF,'(A)') '# Extinction correction'
            #CALL WRVAL(IUCIF,'_gsas_exptl_extinct_corr_T_min',TEXT(1:10))
            #CALL WRVAL(IUCIF,'_gsas_exptl_extinct_corr_T_max',TEXT(11:20))

            if not oneblock:
                # instrumental profile terms go here
                WriteCIFitem('_pd_proc_ls_profile_function','?')

            WriteCIFitem('\n# STRUCTURE FACTOR TABLE')            
            WriteCIFitem('loop_' + 
                         '\n\t' + refprx + 'index_h' + 
                         '\n\t' + refprx + 'index_k' + 
                         '\n\t' + refprx + 'index_l' + 
                         '\n\t_pd_refln_wavelength_id' +
                         '\n\t_pd_refln_phase_id' + 
                         '\n\t' + refprx + 'status' + 
                         '\n\t' + refprx + 'F_squared_meas' + 
                         '\n\t' + refprx + 'F_squared_sigma' + 
                         '\n\t' + refprx + 'F_squared_calc' + 
                         '\n\t' + refprx + 'phase_calc' + 
                         '\n\t_pd_refln_d_spacing' + 
                         '\n\t_gsas_i100_meas')

            WriteCIFitem('_reflns_number_total', text)
            WriteCIFitem('_reflns_limit_h_min', text)
            WriteCIFitem('_reflns_limit_h_max', text)
            WriteCIFitem('_reflns_limit_k_min', text)
            WriteCIFitem('_reflns_limit_k_max', text)
            WriteCIFitem('_reflns_limit_l_min', text)
            WriteCIFitem('_reflns_limit_l_max', text)
            WriteCIFitem('_reflns_d_resolution_high', text)
            WriteCIFitem('_reflns_d_resolution_low', text)
            
            WriteCIFitem('\n# POWDER DATA TABLE')
            # is data fixed step?
            fixedstep = False
            # are ESDs sqrt(I)
            countsdata = False
            zero = 0.01
            if fixedstep:
                WriteCIFitem('_pd_meas_2theta_range_min', text)
                WriteCIFitem('_pd_meas_2theta_range_max', text)
                WriteCIFitem('_pd_meas_2theta_range_inc', text)
                # zero correct
                if zero != 0.0:
                    WriteCIFitem('_pd_proc_2theta_range_min', text)
                    WriteCIFitem('_pd_proc_2theta_range_max', text)
                    WriteCIFitem('_pd_proc_2theta_range_inc', text)
            WriteCIFitem('loop_' +
                         '\n\t_pd_proc_d_spacing')
                         #'_pd_meas_time_of_flight'
            if not fixedstep:
                if zero != 0.0:
                    WriteCIFitem('\t_pd_proc_2theta_corrected')
                else:
                    WriteCIFitem('\t_pd_meas_2theta_scan')
            if countsdata:
                WriteCIFitem('\t_pd_meas_counts_total')
            else:
                WriteCIFitem('\t_pd_meas_intensity_total')
            WriteCIFitem('\t_pd_proc_ls_weight')
            WriteCIFitem('\t_pd_proc_intensity_bkg_calc')
            WriteCIFitem('\t_pd_calc_intensity_total')
            if zero != 0.0:
                WriteCIFitem('_pd_proc_number_of_points', text)
            else:
                WriteCIFitem('_pd_meas_number_of_points', text)
                         
        def WriteSingleXtalData(histlbl):
            histblk = self.GroupedParms['HKLF'][histlbl]
            print 'TODO: single xtal here data for',histblk["Instrument Parameters"][0]['InstrName']
            # see wrreflist.for
            refprx = '_refln.' # mm
            refprx = '_refln_' # normal
            
            WriteCIFitem('\n# STRUCTURE FACTOR TABLE')            
            WriteCIFitem('loop_' + 
                         '\n\t' + refprx + 'index_h' + 
                         '\n\t' + refprx + 'index_k' + 
                         '\n\t' + refprx + 'index_l' + 
                         '\n\t' + refprx + 'status' + 
                         '\n\t' + refprx + 'F_squared_meas' + 
                         '\n\t' + refprx + 'F_squared_sigma' + 
                         '\n\t' + refprx + 'F_squared_calc' + 
                         '\n\t' + refprx + 'phase_calc')
            WriteCIFitem('_reflns_number_total', text)
            WriteCIFitem('_reflns_number_observed', text)
            WriteCIFitem('_reflns_limit_h_min', text)
            WriteCIFitem('_reflns_limit_h_max', text)
            WriteCIFitem('_reflns_limit_k_min', text)
            WriteCIFitem('_reflns_limit_k_max', text)
            WriteCIFitem('_reflns_limit_l_min', text)
            WriteCIFitem('_reflns_limit_l_max', text)
            WriteCIFitem('_reflns_d_resolution_high', text)
            WriteCIFitem('_reflns_d_resolution_low', text)

        #============================================================
        # the export process starts here
        self.loadTree()
        # ====== move this later ==================================================
        #for key in self.OverallParms['Covariance']:
        #    print key
        #    print type(self.OverallParms['Covariance'][key])
        sigDict = dict(zip(
            self.OverallParms['Covariance']['varyList'],
            self.OverallParms['Covariance']['sig'])
                       )
        # get restraint info
        #restraintDict = self.OverallParms.get('Restraints',{})
        for i in  self.OverallParms['Constraints']:
            print i
            for j in self.OverallParms['Constraints'][i]:
                print j
        return
        parmDict = {}
        #G2mv.InitVars()
        # do constraint processing
        #try:
        #    groups,parmlist = G2mv.GroupConstraints(constrDict)
        #    G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList)
        #except:
        #    import GSASIIstrIO as G2stIO
        #    print ' *** ERROR - your constraints are internally inconsistent ***'
        #    errmsg, warnmsg = G2stIO.CheckConstraints(varyList,constrDict,fixedList)
        #    print 'Errors',errmsg
        #    if warnmsg: print 'Warnings',warnmsg
        #    raise Exception(' *** CIF creation aborted ***')
        
        # add the uncertainties into the esd dictionary (sigDict)
        #sigDict.update(G2mv.ComputeDepESD(
        #    self.OverallParms['Covariance']['covMatrix'],
        #    self.OverallParms['Covariance']['varyList'],
        #    parmDict))

        return

        self.CIFdate = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%dT%H:%M")
        # count phases, powder and single crystal histograms
        self.nphase = len(self.GroupedParms.get("Phases",[]))
        self.npowder = len(self.GroupedParms.get("PWDR",[]))
        self.nsingle = len(self.GroupedParms.get("HKLF",[]))
        # is there anything to export?
        if self.nphase + self.npowder + self.nsingle == 0: 
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
        if mode != "full" or self.npowder + self.nsingle == 0:
            self.quickmode = True
            oneblock = True
            if self.nphase == 0:
                self.G2frame.ErrorDialog(
                    'No phase present',
                    'Cannot create a coordinates CIF with no phases')
                return
            elif self.nphase > 1: # quick mode: choose one phase
                choices = sorted(self.GroupedParms['Phases'].keys())
                phasenum = G2gd.ItemSelector(choices,self.G2frame)
                if phasenum is None: return
                phasenam = choices[phasenum]
        # will this require a multiblock CIF?
        elif self.nphase > 1:
            oneblock = False
        elif self.npowder + self.nsingle > 1:
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
            for key in self.GroupedParms:
                if key == 'Phases': continue
                for key1 in self.GroupedParms[key]:
                    if key == "PWDR":
                        key2 = "Sample Parameters"
                        d = self.GroupedParms[key][key1][key2]
                    elif key == "HKLF":
                        key2 = "Instrument Parameters"
                        d = self.GroupedParms[key][key1][key2][0]
                    else:
                        raise Exception,"Unexpected histogram type for CIF: "+str(key)
                        
                    lbllist.append(key1)
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
                if invalid > 1: msg = (
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

        #======================================================================
        # Start writing the CIF - single block
        #======================================================================
        if oneblock:
            if phasenam is None: # if not already selected, select the first phase (should be one) 
                phasenam = self.GroupedParms['Phases'].keys()[0]
            #print 'phasenam',phasenam
            phaseblk = self.GroupedParms['Phases'][phasenam] # pointer to current phase info
            # select data, should only be one set, but take whatever comes 1st
            if not self.quickmode:
                for key in self.GroupedParms:
                    if key == 'Phases': continue
                    for key1 in self.GroupedParms[key]:
                        histblk = self.GroupedParms[key][key1]
                        histkey = key1
                        histtyp = key
                        if key == "PWDR":
                            instnam = histblk["Sample Parameters"]['InstrName']
                        elif key == "HKLF":
                            instnam = histblk["Instrument Parameters"][0]['InstrName']
                        break
                    break
                instnam = instnam.replace(' ','')
            WriteCIFitem('data_'+self.CIFname)
            if not self.quickmode: # publication information
                WriteCIFitem('_pd_block_id',
                             str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                             str(self.shortauthorname) + "|" + instnam)
                WriteAudit()
                WritePubTemplate()
                WriteOverall()
            # report the phase
            WritePhaseTemplate()
            WritePhaseInfo(phasenam)
            if not self.quickmode:
                if histtyp == "PWDR":
                    WritePowderTemplate()
                    WritePowderData(histkey)
                else:
                    WriteSnglXtalTemplate()
                    WriteSingleXtalData(histkey)
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
            if self.nphase > 1:
                loopprefix = ''
                WriteCIFitem('loop_   _pd_phase_block_id')
            else:
                loopprefix = '_pd_phase_block_id'
            
            for i,phasenam in enumerate(sorted(self.GroupedParms['Phases'].keys())):
                datablockidDict[phasenam] = (str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                             'phase_'+ str(i+1) + '|' + str(self.shortauthorname))
                WriteCIFitem(loopprefix,datablockidDict[phasenam])
            # loop over data blocks
            i = 0
            if self.npowder + self.nsingle > 1:
                loopprefix = ''
                WriteCIFitem('loop_   _pd_block_diffractogram_id')
            else:
                loopprefix = '_pd_block_diffractogram_id'
            for key in self.GroupedParms:
                if key == 'Phases': continue
                for key1 in self.GroupedParms[key]:
                    i += 1
                    histblk = self.GroupedParms[key][key1]
                    if key == "PWDR":
                        instnam = histblk["Sample Parameters"]['InstrName']
                    elif key == "HKLF":
                        instnam = histblk["Instrument Parameters"][0]['InstrName']
                    instnam = instnam.replace(' ','')
                    if datablockidDict.get(key1):
                        print "Error: histogram name is duplicated: ",key1
                        return
                    datablockidDict[key1] = (str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                                             str(self.shortauthorname) + "|" +
                                             instnam + "_hist_"+str(i))
                    WriteCIFitem(loopprefix,datablockidDict[key1])
            #============================================================
            # loop over phases, exporting them
            phasebyhistDict = {} # create a cross-reference to phases by histogram
            for j,phasenam in enumerate(sorted(self.GroupedParms['Phases'].keys())):
                i = j + 1
                WriteCIFitem('\ndata_'+self.CIFname+"_phase_"+str(i))
                print "debug, processing ",phasenam
                WriteCIFitem('# Information for phase '+str(i))
                WriteCIFitem('_pd_block_id',datablockidDict[phasenam])
                # report the phase
                WritePhaseTemplate()
                WritePhaseInfo(phasenam)

            #============================================================
            # loop over histograms, exporting them
            i = 0
            for key in self.GroupedParms:
                if key == 'Phases': continue
                for key1 in self.GroupedParms[key]:
                    i += 1
                    histblk = self.GroupedParms[key][key1]
                    if key == "PWDR":
                        WriteCIFitem('\ndata_'+self.CIFname+"_pwd_"+str(i))
                    elif key == "HKLF":
                        WriteCIFitem('\ndata_'+self.CIFname+"_sx_"+str(i))
                    WriteCIFitem('# Information for histogram '+str(i)+': '+
                                 key1)
                    WriteCIFitem('_pd_block_id',datablockidDict[key1])
                    if key == "PWDR":
                        WritePowderTemplate()
                        WritePowderData(key1)
                    else:
                        WriteSnglXtalTemplate()
                        WriteSingleXtalData(key1)

        # TODO: how to report _pd_proc_ls_peak_cutoff?
        WriteCIFitem('#--' + 15*'eof--' + '#')

