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
                if "\n" in value or len(value)+len(name)+4 > 70:
                    if name.strip(): print name
                    print '; '+value
                    print '; '
                elif " " in value:
                    print name,'  ','"' + str(value) + '"'
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
                WriteCIFitem('# Overall powder R-factors')
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
            print 'TODO: phase info for',phasenam,'goes here'
            phasedict = self.GroupedParms['Phases'][phasenam] # pointer to current phase info
            print phasedict.keys()
            # see WRITEPHASE

        def WritePowderData(histlbl):
            histblk = self.GroupedParms['PWDR'][histlbl]
            print 'TODO: powder here data for',histblk["Sample Parameters"]['InstrName']
            # see WRPOWDHIST & WRREFLIST

        def WriteSingleXtalData(histlbl):
            histblk = self.GroupedParms['HKLF'][histlbl]
            print 'TODO: single xtal here data for',histblk["Instrument Parameters"][0]['InstrName']
            # see WRREFLIST

        # the export process starts here
        self.loadTree()
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
            WriteCIFitem('# pointers to phase and histogram blocks')
            # loop over future phase blocks
            if self.nphase > 1:
                loopprefix = ''
                WriteCIFitem('loop_   _pd_phase_block_id')
            else:
                loopprefix = '_pd_phase_block_id'
            for i,phasenam in enumerate(sorted(self.GroupedParms['Phases'].keys())):
                WriteCIFitem(loopprefix,
                             str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                             'phase_'+ str(i) + '|' + str(self.shortauthorname))
            # loop over future data blocks, create and save block names
            i = 0
            datablockidDict = {}
            if self.npowder + self.nsingle > 1:
                loopprefix = ''
                WriteCIFitem('loop_   _pd_block_diffractogramphase_id')
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
                    datablockidDict[key1] = (str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                                             str(self.shortauthorname) + "|" +
                                             instnam + "_hist_"+str(i))
                    WriteCIFitem(loopprefix,datablockidDict[key1])
            #============================================================
            # export phase information
            for i,phasenam in enumerate(sorted(self.GroupedParms['Phases'].keys())):
                WriteCIFitem('\ndata_'+self.CIFname+"_phase_"+str(i))
                print "debug, processing ",phasenam
                WriteCIFitem('# Information for phase '+str(i))
                WriteCIFitem('_pd_block_id',
                             str(self.CIFdate) + "|" + str(self.CIFname) + "|" +
                             'phase_'+ str(i) + '|' + str(self.shortauthorname))
                # pointers to histograms used in this phase
                histlist = []
                for hist in self.GroupedParms['Phases'][phasenam]['Histograms']:
                    if self.GroupedParms['Phases'][phasenam]['Histograms'][hist]['Use']:
                        blockid = datablockidDict.get(hist)
                        if not blockid:
                            print "Internal error: no block for data. Phase "+str(
                                phasenam)+" histogram "+str(hist)
                            histlist = []
                            break
                        histlist.append(blockid)
                if len(histlist) == 0:
                    WriteCIFitem('# Note: phase has no associated data')
                elif len(histlist) == 1:
                    WriteCIFitem('_pd_block_diffractogram_id',histlist[0])
                else:
                    WriteCIFitem('loop_  _pd_block_diffractogram_id')
                    for hist in histlist:
                        WriteCIFitem('',hist)
                # report the phase
                WritePhaseTemplate()
                WritePhaseInfo(phasenam)

            #============================================================
            # loop over histograms
            for key in self.GroupedParms:
                if key == 'Phases': continue
                for key1 in self.GroupedParms[key]:
                    i += 1
                    histblk = self.GroupedParms[key][key1]
                    if key == "PWDR":
                        WriteCIFitem('\ndata_'+self.CIFname+"_p_"+str(i))
                    elif key == "HKLF":
                        WriteCIFitem('\ndata_'+self.CIFname+"_s_"+str(i))
                    WriteCIFitem('# Information for histogram '+str(i)+': '+
                                 key1)
                    WriteCIFitem('_pd_block_id',datablockidDict[key1])
                    if key == "PWDR":
                        WritePowderTemplate()
                        WritePowderData(key1)
                    else:
                        WriteSnglXtalTemplate()
                        WriteSingleXtalData(key1)
        WriteCIFitem('#--' + 15*'eof--' + '#')
