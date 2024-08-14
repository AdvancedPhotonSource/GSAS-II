#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This module initially written by Conrad Gillard. For any enquiries please contact conrad.gillard@gmail.com
# Export3col exporter adapted from Exportbracket by BHT
from __future__ import division, print_function
import wx
import GSASIIpath
import GSASIIIO as G2IO
from collections import OrderedDict
from GSASIImath import ValEsd

class Exportbracket(G2IO.ExportBaseclass):
    '''Enables export of parameters that are commonly needed for publications, in bracket notation
    '''

    def __init__(self, G2frame):
        G2IO.ExportBaseclass.__init__(self,G2frame=G2frame,formatName='Bracket notation CSV',
            extension='.csv',longFormatName='Export commonly needed parameters')
        self.exporttype = ['project']

    def Exporter(self, event=None):
        
        # Define function to extract parameter and sigma from covariances, for later use
        def GetParamSig(phase_num, hist_num, keyword, display_name):
            param_index = None
            try:
                param_index = self.OverallParms['Covariance']['varyList'].index(phase_num + ':' + hist_num + keyword)
            except:
                pass
            if param_index is not None:
                param = self.OverallParms['Covariance']['variables'][param_index]
                # Extract parameter uncertainty
                param_sig = self.OverallParms['Covariance']['sig'][param_index]
                # Create dictionary entry containing parameter and sigma in bracket notation
                model_parameters[display_name] = ValEsd(param, param_sig)

        # # set up for export
        self.InitExport(event)
        if self.ExportSelect(): return  # set export parameters; get file name
        self.OpenFile()
        wx.BeginBusyCursor()

        # Export model parameters in bracket notation
        try:
            # Initialise ordered dictionary to hold model parameters and uncertainties
            model_parameters = OrderedDict()

            # load all of the tree into a set of dicts
            self.loadTree()
            # create a dict with refined values and their uncertainties
            self.loadParmDict()

            # Initialise phase counter, for later use
            phase_num = 0
            # Extract lattice parameters and uncertainties, and convert to bracket notation
            for phasedict in self.Phases.items():
                phasenam = phasedict[0]
                cellList, cellSig = self.GetCell(phasenam)

                # Set up list of cell quantity symbols and units, to be used for generating text in output CSV
                cell_quantity_symbols = ["a", "b", "c", "alpha", "beta", "gamma", "Volume"]
                cell_quantity_units = ["(Å)", "(Å)", "(Å)", "(°)", "(°)", "(°)", "(Å³)"]

                for i in range(0, len(cellList)):
                    if cellSig[i] > 0:
                        # Formulate lattice parameter in bracket notation
                        current_lp_bracket = ValEsd(cellList[i], cellSig[i])
                        # Write to dictionary that will later be exported to CSV
                        model_parameters[phasenam + " " + cell_quantity_symbols[i] + " " +
                                         cell_quantity_units[i]] = current_lp_bracket

                # Get phase and weight fractions uncertainties, if they have been refined
                for hist_num,hist_name in enumerate(self.Histograms):
                    try:
                        # Get phase fraction uncertainty, if phase fractions have been refined
                        phasefrac_unc = self.sigDict[str(phase_num) + ':' + str(hist_num) + ':Scale']
                        # Get name of histogram associated with this phase, for later use
                        #hist_name = list(self.Histograms.keys())[hist_num]
                        # Extract phase fraction value
                        phasefrac = phasedict[1]['Histograms'][hist_name]['Scale'][0]
                        # Write phase if there is more than one histogram, specify which one
                        if len(self.Histograms) > 1:
                            model_parameters[phasenam + " Phase Fraction in: " + hist_name] = (
                                ValEsd(phasefrac, phasefrac_unc))
                        # If there is only one histogram, no need to specify which histogram the fraction is based on
                        else:
                            model_parameters[phasenam + " Phase Fraction"] = ValEsd(phasefrac, phasefrac_unc)
                    except:
                        pass

                    try:
                        var = str(phase_num) + ':' + str(hist_num) + ':WgtFrac'
                        depSigDict = self.OverallParms['Covariance'].get('depSigDict',{})
                        weight_frac,weight_frac_unc = depSigDict.get(var,[0,None])

                        # Write phase + weight fractions in bracket notation to dictionary, to be exported as a CSV
                        # If there is more than one histogram, specify which one the fraction is based on
                        if len(self.Histograms) > 1:
                            model_parameters[phasenam + " Weight Fraction in: " + hist_name] = (
                                ValEsd(weight_frac, weight_frac_unc))
                        # If there is only one histogram, no need to specify which histogram the fraction is based on
                        else:
                            model_parameters[phasenam + " Weight Fraction"] = ValEsd(weight_frac, weight_frac_unc)
                    except:
                        pass

                    # Get preferred orientation details for phase, if refined
                    try:
                        pref_orr_props = phasedict[1]['Histograms'][hist_name]['Pref.Ori.']
                        # Check if March Dollase has been refined
                        if pref_orr_props[2] and pref_orr_props[0] == "MD":
                           # If so, first extract MD axis and write to dictionary to be exported
                           MD_axis = "\'" + ''.join(map(str, pref_orr_props[3]))
                           model_parameters[phasenam + " March Dollase Axis"] = MD_axis
                           # Next, extract MD ratio
                           MD_ratio = pref_orr_props[1]
                           MD_sig = self.sigDict[str(phase_num)+":" + str(hist_num) + ":MD"]
                           # Formulate MD ratio in bracket notation
                           MD_bracket = ValEsd(MD_ratio, MD_sig)
                           # Write MD ratio to dictionary to be exported
                           model_parameters[phasenam + " March Dollase Ratio"] = MD_bracket
                    except:
                        pass
                # Increment phase number counter
                phase_num += 1

# Extract sample displacements, zero offset and D(ij)s (if refined)
            for i,hist_name in enumerate(self.Histograms):
                hist_num = str(i)
                # Extract zero offset, if refined
                GetParamSig("", hist_num, ':Zero', "Zero Offset")
                # Extract Bragg-Brentano sample displacement, if refined
                GetParamSig("", hist_num, ':Shift', "Sample Displacement (micron)")
                # Extract Debye-Scherrer sample X displacement, if refined
                GetParamSig("", hist_num, ':DisplaceX', "Sample X Displacement (micron)")
                # Extract Debye-Scherrer sample Y displacement, if refined
                GetParamSig("", hist_num, ':DisplaceY', "Sample Y Displacement (micron)")
                # Extract hydrostatic strains, if refined
                for phase_num,phase_name in enumerate(self.Phases):
                    for d_i in range(1, 4):
                        for d_j in range(1, 4):
                            GetParamSig(str(phase_num), hist_num, ':D' + str(d_i) + str(d_j),
                                        phase_name + ' D' + str(d_i) + str(d_j))

                # Extract atomic parameters, if refined
                for phase_num,phase_name in enumerate(self.Phases):
                    # atom_list = list(self.Phases.values())[phase_num]["Atoms"]
                    atom_list = self.Phases[phase_name]["Atoms"] #  same as above?
                    for atom_num,atom in enumerate(atom_list):
                        # Extract isotropic thermal parameters, if refined
                        GetParamSig(str(phase_num), ':', 'AUiso:' + str(atom_num),
                                    phase_name + ' ' + atom[0] + ' Uiso')
                        # Extract anisotropic thermal parameters (Uijs), if refined
                        for Ui in range(1, 4):
                            for Uj in range(1, 4):
                                GetParamSig(str(phase_num), ':', 'AU' + str(Ui) + str(Uj) + ':' + str(atom_num),
                                            phase_name + ' ' + atom[0] + ' U' + str(Ui) + str(Uj))
                        # Extract fractional occupancies, if refined
                        GetParamSig(str(phase_num), ':', 'Afrac:' + str(atom_num),
                                    phase_name + ' ' + atom[0] + ' Occupancy')
                        # Extract atom X Y Z, if refined
                        for atom_axis in ('x', 'y', 'z'):
                            variable_code = str(phase_num) + ':' + ':' + 'dA' + atom_axis + ':' + str(atom_num)
                            # Get sigma, if refined
                            try:
                                param_index = self.OverallParms['Covariance']['varyList'].index(variable_code)
                                # Extract uncertainty
                                atom_axis_sig = self.OverallParms['Covariance']['sig'][param_index]
                                # Extract value
                                atom_axis_val = list(self.Phases.values())[phase_num]["Atoms"][atom_num][ord(atom_axis)-117]
                                # Convert to bracket notation and add to dictionary, which will be exported as a CSV
                                model_parameters[phase_name + ' ' + atom[0] + ' ' + atom_axis] = \
                                    ValEsd(atom_axis_val, atom_axis_sig)
                            except: pass

            # Extract rWp
            rWP = self.OverallParms['Covariance']['Rvals']['Rwp']
            # Write to dictionary to be printed, rounding to 3 significant figures for readability
            model_parameters["wR"] = str(ValEsd(rWP, -0.1)) + "%"

            # Write to CSV
            # parameter_names = ""
            # for parameter_name in model_parameters.keys():
            #     parameter_names = parameter_names + str(parameter_name) + ", "
            # self.Write(parameter_names[0:-2])

            # parameter_values = ""
            # for parameter_value in model_parameters.values():
            #     parameter_values = parameter_values + str(parameter_value) + ", "
            # self.Write(parameter_values[0:-2])
            
            for name in model_parameters:
                self.Write('%s, %s,'%(name,model_parameters[name]))

        finally:
            wx.EndBusyCursor()
        self.CloseFile()


class Export3col(G2IO.ExportBaseclass):
    '''Enables export of parameters that are commonly needed for publications, with esds
    in a separate column
    '''

    def __init__(self, G2frame):
        G2IO.ExportBaseclass.__init__(self,G2frame=G2frame,formatName='common prm CSV',
            extension='.csv',longFormatName='Export commonly needed parameters with s.u. in a separate column')
        self.exporttype = ['project']

    def ValEsd2col(self, param, param_sig):
        '''Return two values with the formated value as the first number and the 
        standard uncertainty (if provided) as the second value.
        '''
        col1 = ValEsd(param, -abs(param_sig))
        col2 = ''
        if param_sig > 0:
            col2 = ValEsd(param_sig, -param_sig/100)
        return col1,col2            
        
    def Exporter(self, event=None):
        
        # Define function to extract parameter and sigma from covariances, for later use
        def GetParamSig(phase_num, hist_num, keyword, display_name):
            param_index = None
            try:
                param_index = self.OverallParms['Covariance']['varyList'].index(phase_num + ':' + hist_num + keyword)
            except:
                pass
            if param_index is not None:
                param = self.OverallParms['Covariance']['variables'][param_index]
                # Extract parameter uncertainty
                param_sig = self.OverallParms['Covariance']['sig'][param_index]
                # Create dictionary entry containing parameter and sigma in bracket notation
                model_parameters[display_name] = self.ValEsd2col(param, param_sig)

        # # set up for export
        self.InitExport(event)
        if self.ExportSelect(): return  # set export parameters; get file name
        self.OpenFile()
        wx.BeginBusyCursor()

        # Export model parameters in bracket notation
        try:
            # Initialise ordered dictionary to hold model parameters and uncertainties
            model_parameters = OrderedDict()

            # load all of the tree into a set of dicts
            self.loadTree()
            # create a dict with refined values and their uncertainties
            self.loadParmDict()

            # Initialise phase counter, for later use
            phase_num = 0
            # Extract lattice parameters and uncertainties, and convert to bracket notation
            for phasedict in self.Phases.items():
                phasenam = phasedict[0]
                cellList, cellSig = self.GetCell(phasenam)
                # Initialise lattice parameter letter
                lp_letter = "a"
                for i in range(0, len(cellList)):
                # for cell in cellList:
                    if cellSig[i] > 0:
                        # Formulate lattice parameter 
                        model_parameters[phasenam + " " + lp_letter + " (Å)"] = self.ValEsd2col(cellList[i], cellSig[i])
                        # Increment lattice parameter letter
                        lp_letter = chr(ord(lp_letter[0]) + 1)
                    else:
                        break

                # Get phase and weight fractions uncertainties, if they have been refined
                for hist_num,hist_name in enumerate(self.Histograms):
                    try:
                        # Get phase fraction uncertainty, if phase fractions have been refined
                        phasefrac_unc = self.sigDict[str(phase_num) + ':' + str(hist_num) + ':Scale']
                        # Get name of histogram associated with this phase, for later use
                        #hist_name = list(self.Histograms.keys())[hist_num]
                        # Extract phase fraction value
                        phasefrac = phasedict[1]['Histograms'][hist_name]['Scale'][0]
                        # Write phase if there is more than one histogram, specify which one
                        if len(self.Histograms) > 1:
                            model_parameters[phasenam + " Phase Fraction in: " + hist_name] = (
                                self.ValEsd2col(phasefrac, phasefrac_unc))
                        # If there is only one histogram, no need to specify which histogram the fraction is based on
                        else:
                            model_parameters[phasenam + " Phase Fraction"] = self.ValEsd2col(phasefrac, phasefrac_unc)
                    except:
                        pass

                    try:
                        var = str(phase_num) + ':' + str(hist_num) + ':WgtFrac'
                        depSigDict = self.OverallParms['Covariance'].get('depSigDict',{})
                        weight_frac,weight_frac_unc = depSigDict.get(var,[0,None])

                        # Write phase + weight fractions in bracket notation to dictionary, to be exported as a CSV
                        # If there is more than one histogram, specify which one the fraction is based on
                        if len(self.Histograms) > 1:
                            model_parameters[phasenam + " Weight Fraction in: " + hist_name] = (
                                self.ValEsd2col(weight_frac, weight_frac_unc))
                        # If there is only one histogram, no need to specify which histogram the fraction is based on
                        else:
                            model_parameters[phasenam + " Weight Fraction"] = self.ValEsd2col(weight_frac, weight_frac_unc)
                    except:
                        pass

                    # Get preferred orientation details for phase, if refined
                    try:
                        pref_orr_props = phasedict[1]['Histograms'][hist_name]['Pref.Ori.']
                        # Check if March Dollase has been refined
                        if pref_orr_props[2] and pref_orr_props[0] == "MD":
                           # If so, first extract MD axis and write to dictionary to be exported
                           MD_axis = "\'" + ''.join(map(str, pref_orr_props[3]))
                           model_parameters[phasenam + " March Dollase Axis"] = (MD_axis,'')
                           # Next, extract MD ratio
                           MD_ratio = pref_orr_props[1]
                           MD_sig = self.sigDict[str(phase_num)+":" + str(hist_num) + ":MD"]
                           # Formulate MD ratio in bracket notation
                           MD_bracket = self.ValEsd2col(MD_ratio, MD_sig)
                           # Write MD ratio to dictionary to be exported
                           model_parameters[phasenam + " March Dollase Ratio"] = (MD_bracket,'')
                    except:
                        pass
                # Increment phase number counter
                phase_num += 1

# Extract sample displacements, zero offset and D(ij)s (if refined)
            for i,hist_name in enumerate(self.Histograms):
                hist_num = str(i)
                # Extract zero offset, if refined
                GetParamSig("", hist_num, ':Zero', "Zero Offset")
                # Extract Bragg-Brentano sample displacement, if refined
                GetParamSig("", hist_num, ':Shift', "Sample Displacement (micron)")
                # Extract Debye-Scherrer sample X displacement, if refined
                GetParamSig("", hist_num, ':DisplaceX', "Sample X Displacement (micron)")
                # Extract Debye-Scherrer sample Y displacement, if refined
                GetParamSig("", hist_num, ':DisplaceY', "Sample Y Displacement (micron)")
                # Extract hydrostatic strains, if refined
                for phase_num,phase_name in enumerate(self.Phases):
                    for d_i in range(1, 4):
                        for d_j in range(1, 4):
                            GetParamSig(str(phase_num), hist_num, ':D' + str(d_i) + str(d_j),
                                        phase_name + ' D' + str(d_i) + str(d_j))

                # Extract atomic parameters, if refined
                for phase_num,phase_name in enumerate(self.Phases):
                    # atom_list = list(self.Phases.values())[phase_num]["Atoms"]
                    atom_list = self.Phases[phase_name]["Atoms"] #  same as above?
                    for atom_num,atom in enumerate(atom_list):
                        # Extract isotropic thermal parameters, if refined
                        GetParamSig(str(phase_num), ':', 'AUiso:' + str(atom_num),
                                    phase_name + ' ' + atom[0] + ' Uiso')
                        # Extract anisotropic thermal parameters (Uijs), if refined
                        for Ui in range(1, 4):
                            for Uj in range(1, 4):
                                GetParamSig(str(phase_num), ':', 'AU' + str(Ui) + str(Uj) + ':' + str(atom_num),
                                            phase_name + ' ' + atom[0] + ' U' + str(Ui) + str(Uj))
                        # Extract fractional occupancies, if refined
                        GetParamSig(str(phase_num), ':', 'Afrac:' + str(atom_num),
                                    phase_name + ' ' + atom[0] + ' Occupancy')
                        # Extract atom X Y Z, if refined
                        for atom_axis in ('x', 'y', 'z'):
                            variable_code = str(phase_num) + ':' + ':' + 'dA' + atom_axis + ':' + str(atom_num)
                            # Get sigma, if refined
                            try:
                                param_index = self.OverallParms['Covariance']['varyList'].index(variable_code)
                                # Extract uncertainty
                                atom_axis_sig = self.OverallParms['Covariance']['sig'][param_index]
                                # Extract value
                                atom_axis_val = list(self.Phases.values())[phase_num]["Atoms"][atom_num][ord(atom_axis)-117]
                                # Convert to bracket notation and add to dictionary, which will be exported as a CSV
                                model_parameters[phase_name + ' ' + atom[0] + ' ' + atom_axis] = \
                                    self.ValEsd2col(atom_axis_val, atom_axis_sig)
                            except: pass

            # Extract rWp
            rWP = self.OverallParms['Covariance']['Rvals']['Rwp']
            # Write to dictionary to be printed, rounding to 3 significant figures for readability
            model_parameters["wR"] = (ValEsd(rWP, -0.1) + '%', '')

            # Write to CSV
            for name in model_parameters:
                self.Write('{:}, {:}, {:}'.format(name,*model_parameters[name]))

        finally:
            wx.EndBusyCursor()
        self.CloseFile()
        
