# -*- coding: utf-8 -*-
'''NIST XRD Fundamental Parameters interface routines follow:
'''
from __future__ import division, print_function
import os.path
import numpy as np
import copy

import wx
import wx.lib.scrolledpanel as wxscroll

from . import NIST_profile as FP

from . import GSASIIpath
from . import GSASIIctrlGUI as G2G
from . import GSASIIdataGUI as G2gd
from . import GSASIIplot as G2plt
from . import GSASIImath as G2mth
from . import GSASIIpwd as G2pwd
from . import GSASIIfiles as G2fil
WACV = wx.ALIGN_CENTER_VERTICAL

simParms = {}
'''Parameters to set range for pattern simulation
'''

parmDict = {}
'''Parameter dict used for reading Topas-style values. These are
converted to SI units and placed into :data:`NISTparms`
'''

NISTparms = {}
'''Parameters in a nested dict, with an entry for each concolutor. Entries in
those dicts have values in SI units (of course). NISTparms can be
can be input directly or can be from created from :data:`parmDict`
by :func:`XferFPAsettings`
'''

BraggBrentanoParms = [
    ('divergence', 0.5, 'Bragg-Brentano divergence angle (degrees)'),
    ('soller_angle', 2.0, 'Soller slit axial divergence (degrees)'),
    ('Rs', 220, 'Diffractometer radius (mm)'),
    ('filament_length', 12., 'X-ray tube line focus length (mm)'),
    ('sample_length', 12., 'Illuminated sample length in axial direction (mm)'),
    ('receiving_slit_length', 12., 'Length of receiving slit in axial direction (mm)'),
    ('LAC_cm', 0.,'Linear absorption coef. adjusted for packing density (cm-1)'),
    ('sample_thickness', 1., 'Depth of sample (mm)'),
    ('convolution_steps', 8, 'Number of Fourier-space bins per two-theta step'),
    ('source_width', 0.04,'Tube filament width, in projection at takeoff angle (mm)'),
    ('tube-tails_L-tail', -1.,'Left-side tube tails width, in projection (mm)'),
    ('tube-tails_R-tail', 1.,'Right-side tube tails width, in projection (mm)'),
    ('tube-tails_rel-I', 0.001,'Tube tails fractional intensity (no units)'),
    ]
'''FPA dict entries used in :func:`FillParmSizer`. Tuple contains
a dict key, a default value and a description. These are the parameters
needed for all Bragg Brentano instruments
'''

BBPointDetector = [
    ('receiving_slit_width', 0.2, 'Width of receiving slit (mm)'),]
'''Additional FPA dict entries used in :func:`FillParmSizer`
needed for Bragg Brentano instruments with point detectors.
'''

BBPSDDetector = [
    ('SiPSD_th2_angular_range', 3.0, 'Angular range observed by PSD (degrees 2Theta)'),]
'''Additional FPA dict entries used in :func:`FillParmSizer`
needed for Bragg Brentano instruments with linear (1-D) Si PSD detectors.
'''

IBmonoParms = [
    ('src_mono_mm',119,'Distance from xray line source to monochromator crystal (mm)'),
    ('focus_mono_mm',217,'Distance from monochromator crystal to focus slit (mm)'),
    ('passband_mistune',-0.145,'Offset from the tuning of the IBM to the center of the reference line of the spectrum, in units of its bandwidth'),
    ('mono_src_proj_mn',51,'Projection width of line-focus xray source on the monochromator along beam direction (microns), sets bandwidth'),
    ('passband_shoulder',.087,'Width of transition region from flat top to tails, in units of the bandwidth'),
    ('two_theta_mono',27.27,'The full diffraction angle of the IBM crystal, e.g. double 2theta-Bragg for the mono (deg)'),
    ('mono_slit_attenuation',0.03,'Attenuation of Cu K alpha2 lines due to focal slit'),
    ]
'''Additional FPA dict entries used in :func:`FillParmSizer`, needed for Incident Beam Monochromator
'''

IBmono = False
'''set to True if an incident beam monochromator is in use
'''
DetMode = 'BBpoint'
'''The type of detector, either 'BBpoint' for Bragg-Brentano point detector or
      or 'BBPSD' (linear) position sensitive detector
'''

def SetCu2Wave():
    '''Set the parameters to the two-line Cu K alpha 1+2 spectrum
    '''
    parmDict['wave'] = {i:v for i,v in enumerate((1.5405925, 1.5443873))}
    parmDict['int'] = {i:v for i,v in enumerate((0.653817, 0.346183))}
    parmDict['lwidth'] = {i:v for i,v in enumerate((0.501844,0.626579))}

def SetCu6wave():
    '''Set the emission parameters to the NIST six-line Cu K alpha spectrum
    '''
    # values from Marcus Mendenhall from atan_windowed_FP_profile.py
    parmDict['wave'] = {i:v for i,v in enumerate((1.5405925, 1.5443873, 1.5446782, 1.5410769, 1.53471, 1.53382, ))}
    parmDict['int'] = {i:v for i,v in enumerate((0.58384351, 0.2284605 , 0.11258773, 0.07077796, 0.0043303, 0.00208613, ))}
    parmDict['lwidth'] = {i:v for i,v in enumerate((0.436, 0.487, 0.63, 0.558, 2.93, 2.93,))}

def SetMonoWave():
    '''Eliminates the short-wavelength line from the six-line Cu K
    alpha spectrum when incident beam mono; resets it to 6 if no mono
    '''
    if IBmono and len(parmDict['wave']) == 6:
        for key in 'wave','int','lwidth':
            if 5 in parmDict[key]: del parmDict[key][5]
            if 4 in parmDict[key]: del parmDict[key][4]
    if (not IBmono) and len(parmDict['wave']) == 4:
        SetCu6wave()

def writeNIST(filename):
    '''Write the NIST FPA terms into a JSON-like file that can be reloaded
    in _onReadFPA
    '''
    if not filename: return

    fp = open(filename,'w')
    fp.write('# parameters to be used in the NIST XRD Fundamental Parameters program\n')
    fp.write('{\n')
    for key in sorted(NISTparms):
        fp.write("   '"+key+"' : "+str(NISTparms[key])+",")
        if not key: fp.write('  # global parameters')
        fp.write('\n')
    fp.write('}\n')
    fp.close()

#SetCu2Wave() # use these as default
SetCu6wave() # use these as default
SetMonoWave()

def FillParmSizer():
    '''Create a list of input items for the parameter section of the
    input window, sets default values when not set and displays them
    in the scrolledpanel prmPnl.
    '''
    prmSizer = prmPnl.GetSizer()
    prmSizer.Clear(True)
    if IBmono:
        itemList = [i for i in BraggBrentanoParms if not i[0].startswith('tube-tails')]
    else:
        itemList = copy.deepcopy(BraggBrentanoParms)
    if DetMode == 'BBpoint':
        itemList += BBPointDetector
    elif DetMode == 'BBPSD':
        itemList += BBPSDDetector
    else:
        raise Exception('Unknown DetMode in FillParmSizer: '+DetMode)
    if IBmono:
        itemList += IBmonoParms
    text = wx.StaticText(prmPnl,wx.ID_ANY,'label',style=wx.ALIGN_CENTER,
                size=(170,-1)) # make just a bit bigger than largest item in column
    text.SetBackgroundColour(wx.WHITE)
    text.SetForegroundColour(wx.BLACK)
    prmSizer.Add(text,0,wx.EXPAND)
    text = wx.StaticText(prmPnl,wx.ID_ANY,'value',style=wx.ALIGN_CENTER)
    text.SetBackgroundColour(wx.WHITE)
    text.SetForegroundColour(wx.BLACK)
    prmSizer.Add(text,0,wx.EXPAND)
    txtSizer = wx.BoxSizer(wx.HORIZONTAL)
    text = wx.StaticText(prmPnl,wx.ID_ANY,'explanation',style=wx.ALIGN_CENTER)
    text.SetBackgroundColour(wx.WHITE)
    text.SetForegroundColour(wx.BLACK)
    txtSizer.Add(text,1,wx.EXPAND)
    txtSizer.Add(G2G.HelpButton(prmPnl,helpIndex='FPAinput'))
    prmSizer.Add(txtSizer,0,wx.EXPAND)
    for lbl,defVal,text in itemList:
        prmSizer.Add(wx.StaticText(prmPnl,wx.ID_ANY,lbl),1,wx.ALIGN_RIGHT|WACV,1)
        if lbl not in parmDict: parmDict[lbl] = defVal
        ctrl = G2G.ValidatedTxtCtrl(prmPnl,parmDict,lbl,size=(70,-1))
        prmSizer.Add(ctrl,1,wx.ALL|WACV,1)
        txt = wx.StaticText(prmPnl,wx.ID_ANY,text,size=(400,-1))
        txt.Wrap(380)
        prmSizer.Add(txt)
    px1,py = prmSizer.GetSize()
    prmSizer.Layout()
    FPdlg = prmPnl.GetParent()
    FPdlg.SendSizeEvent()

def MakeTopasFPASizer(G2frame,FPdlg,SetButtonStatus):
    '''Create a GUI with parameters for the NIST XRD Fundamental Parameters Code.
    Parameter input is modeled after Topas input parameters.

    :param wx.Frame G2frame: main GSAS-II window
    :param wx.Window FPdlg: Frame or Dialog where GUI will appear
    :param SetButtonStatus: a callback function to call to see what buttons in
      this windows can be enabled. Called with done=True to trigger closing
      the parent window as well.
    :returns: a sizer with the GUI controls
    '''
    def _onOK(event):
        XferFPAsettings(parmDict)
        SetButtonStatus(done=True) # done=True triggers the simulation
        FPdlg.Destroy()
    def _onClose(event):
        SetButtonStatus()
        FPdlg.Destroy()
    def _onAddWave(event):
        newkey = max(parmDict['wave'].keys())+1
        for key,defVal in zip(
            ('wave','int','lwidth'),
            (0.0,   1.0,   0.1),
            ):
            parmDict[key][newkey] = defVal
        wx.CallAfter(MakeTopasFPASizer,G2frame,FPdlg,SetButtonStatus)
    def _onRemWave(event):
        lastkey = max(parmDict['wave'].keys())
        for key in ('wave','int','lwidth'):
            if lastkey in parmDict[key]:
                del parmDict[key][lastkey]
        wx.CallAfter(MakeTopasFPASizer,G2frame,FPdlg,SetButtonStatus)
    def _onSetCu6wave(event):
        SetCu6wave()
        SetMonoWave()
        wx.CallAfter(MakeTopasFPASizer,G2frame,FPdlg,SetButtonStatus)
    def _onSetCu2Wave(event):
        SetCu2Wave()
        wx.CallAfter(MakeTopasFPASizer,G2frame,FPdlg,SetButtonStatus)
    def _onSetDetBtn(event):
        global DetMode
        if detBtn1.GetValue():
            DetMode = 'BBpoint'
            wx.CallAfter(FillParmSizer)
        else:
            DetMode = 'BBPSD'
            wx.CallAfter(FillParmSizer)
    def _onSetMonoBtn(event):
        global IBmono
        IBmono = not monoBtn1.GetValue()
        SetMonoWave()
        #wx.CallAfter(FillParmSizer)
        wx.CallAfter(MakeTopasFPASizer,G2frame,FPdlg,SetButtonStatus)
    def PlotTopasFPA(event):
        XferFPAsettings(parmDict)
        ttArr = np.arange(max(0.5,
                            simParms['plotpos']-simParms['calcwid']),
                            simParms['plotpos']+simParms['calcwid'],
                            simParms['step'])
        intArr = np.zeros_like(ttArr)
        NISTpk = setupFPAcalc()
        try:
            center_bin_idx,peakObj = doFPAcalc(
                NISTpk,ttArr,simParms['plotpos'],simParms['calcwid'],
                simParms['step'])
        except Exception as err:
            msg = "Error computing convolution, revise input"
            print(msg)
            print(err)
            return
        G2plt.PlotFPAconvolutors(G2frame,NISTpk)
        pkPts = len(peakObj.peak)
        pkMax = peakObj.peak.max()
        startInd = center_bin_idx-(pkPts//2) #this should be the aligned start of the new data
        # scale peak so max I=10,000 and add into intensity array
        if startInd < 0:
            intArr[:startInd+pkPts] += 10000 * peakObj.peak[-startInd:]/pkMax
        elif startInd > len(intArr):
            return
        elif startInd+pkPts >= len(intArr):
            offset = pkPts - len( intArr[startInd:] )
            intArr[startInd:startInd+pkPts-offset] += 10000 * peakObj.peak[:-offset]/pkMax
        else:
            intArr[startInd:startInd+pkPts] += 10000 * peakObj.peak/pkMax
        G2plt.PlotXY(G2frame, [(ttArr, intArr)],
                     labelX=r'$2\theta, deg$',
                     labelY=r'Intensity (arbitrary)',
                     Title='FPA peak', newPlot=True, lines=True)
    def _onSaveFPA(event):
        XferFPAsettings(parmDict)
        filename = G2G.askSaveFile(G2frame,'','.NISTfpa',
                                       'dict of NIST FPA values',FPdlg)
        writeNIST(filename)

    if FPdlg.GetSizer(): FPdlg.GetSizer().Clear(True)
    MainSizer = wx.BoxSizer(wx.VERTICAL)
    MainSizer.Add((-1,5))
    waveSizer = wx.FlexGridSizer(cols=len(parmDict['wave'])+1,hgap=3,vgap=5)
    for lbl,prm,defVal in zip(
            (u'Wavelength (\u212b)','Rel. Intensity',u'Lorentz Width\n(\u212b/1000)'),
            ('wave','int','lwidth'),
            (0.0,   1.0,   0.1),
            ):
        text = wx.StaticText(FPdlg,wx.ID_ANY,lbl,style=wx.ALIGN_CENTER)
        text.SetBackgroundColour(wx.WHITE)
        text.SetForegroundColour(wx.BLACK)
        waveSizer.Add(text,0,wx.EXPAND)
        if prm not in parmDict: parmDict[prm] = {}
        for i in parmDict['wave'].keys():
            if i not in parmDict[prm]: parmDict[prm][i] = defVal
            if prm == 'wave':
                ctrl = G2G.ValidatedTxtCtrl(FPdlg,parmDict[prm],i,size=(90,-1),nDig=(10,6))
            else:
                ctrl = G2G.ValidatedTxtCtrl(FPdlg,parmDict[prm],i,size=(90,-1),nDig=(10,3))
            waveSizer.Add(ctrl,1,WACV,1)
    MainSizer.Add(waveSizer)
    MainSizer.Add((-1,5))
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = wx.Button(FPdlg, wx.ID_ANY,'Add wave')
    btnsizer.Add(btn)
    btn.Bind(wx.EVT_BUTTON,_onAddWave)
    btn = wx.Button(FPdlg, wx.ID_ANY,'Remove wave')
    btnsizer.Add(btn)
    btn.Bind(wx.EVT_BUTTON,_onRemWave)
    btn = wx.Button(FPdlg, wx.ID_ANY,'CuKa1+2')
    btnsizer.Add(btn)
    btn.Bind(wx.EVT_BUTTON,_onSetCu2Wave)
    btn = wx.Button(FPdlg, wx.ID_ANY,'NIST CuKa')
    btnsizer.Add(btn)
    btn.Bind(wx.EVT_BUTTON,_onSetCu6wave)
    MainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER, 0)
    MainSizer.Add((-1,5))

    btnsizer = wx.GridBagSizer( 2, 5)
    btnsizer.Add( wx.StaticText(FPdlg, wx.ID_ANY, 'Detector type'),
                 (0,0), (2,1), wx.ALIGN_CENTER | wx.ALL, 5)
    detBtn1 = wx.RadioButton(FPdlg,wx.ID_ANY,'Point',style=wx.RB_GROUP)
    detBtn1.SetValue(DetMode == 'BBpoint')
    btnsizer.Add(detBtn1, (0,1))
    detBtn1.Bind(wx.EVT_RADIOBUTTON,_onSetDetBtn)
    detBtn2 = wx.RadioButton(FPdlg,wx.ID_ANY,'PSD')
    detBtn2.SetValue(not DetMode == 'BBpoint')
    btnsizer.Add(detBtn2, (1,1))
    detBtn2.Bind(wx.EVT_RADIOBUTTON,_onSetDetBtn)
    btnsizer.Add( (40,-1), (0,2), (1,1), wx.ALIGN_CENTER | wx.ALL, 5)
    btnsizer.Add( wx.StaticText(FPdlg, wx.ID_ANY, 'Incident Beam Mono'),
                 (0,3), (2,1), wx.ALIGN_CENTER | wx.ALL, 5)
    monoBtn1 = wx.RadioButton(FPdlg,wx.ID_ANY,'No',style=wx.RB_GROUP)
    monoBtn1.SetValue(not IBmono)
    btnsizer.Add(monoBtn1, (0,4))
    monoBtn1.Bind(wx.EVT_RADIOBUTTON,_onSetMonoBtn)
    monoBtn2 = wx.RadioButton(FPdlg,wx.ID_ANY,'Yes')
    monoBtn2.SetValue(IBmono)
    btnsizer.Add(monoBtn2, (1,4))
    monoBtn2.Bind(wx.EVT_RADIOBUTTON,_onSetMonoBtn)
    MainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER, 0)

    global prmPnl
    prmPnl = wxscroll.ScrolledPanel(FPdlg, wx.ID_ANY, #size=(200,200),
        style = wx.TAB_TRAVERSAL|wx.SUNKEN_BORDER)
    prmSizer = wx.FlexGridSizer(cols=3,hgap=3,vgap=5)
    prmPnl.SetSizer(prmSizer)
    FillParmSizer()
    MainSizer.Add(prmPnl,1,wx.EXPAND,1)
    prmPnl.SetAutoLayout(1)
    prmPnl.SetupScrolling()

    MainSizer.Add((-1,4))
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = wx.Button(FPdlg, wx.ID_ANY, 'Plot peak')
    btnsizer.Add(btn)
    btn.Bind(wx.EVT_BUTTON,PlotTopasFPA)
    btnsizer.Add(wx.StaticText(FPdlg,wx.ID_ANY,' at '))
    if 'plotpos' not in simParms: simParms['plotpos'] =  simParms['minTT']
    ctrl = G2G.ValidatedTxtCtrl(FPdlg,simParms,'plotpos',size=(70,-1))
    btnsizer.Add(ctrl)
    btnsizer.Add(wx.StaticText(FPdlg,wx.ID_ANY,' deg.  '))
    saveBtn = wx.Button(FPdlg, wx.ID_ANY,'Save FPA dict')
    btnsizer.Add(saveBtn)
    saveBtn.Bind(wx.EVT_BUTTON,_onSaveFPA)
    MainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER, 0)
    MainSizer.Add((-1,4))
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    OKbtn = wx.Button(FPdlg, wx.ID_OK)
    OKbtn.SetDefault()
    btnsizer.Add(OKbtn)
    Cbtn = wx.Button(FPdlg, wx.ID_CLOSE,"Cancel")
    btnsizer.Add(Cbtn)
    MainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER, 0)
    MainSizer.Add((-1,4))
    # bindings for close of window
    OKbtn.Bind(wx.EVT_BUTTON,_onOK)
    Cbtn.Bind(wx.EVT_BUTTON,_onClose)
    FPdlg.SetSizer(MainSizer)
    MainSizer.Layout()
    MainSizer.Fit(FPdlg)
    # control window size
    px,py = prmSizer.GetSize()
    dx,dy = FPdlg.GetSize()
    FPdlg.SetMinSize((-1,-1))
    FPdlg.SetMaxSize((-1,-1))
    FPdlg.SetMinSize((dx,dy+200)) # leave a min of 200 points for scroll panel
    FPdlg.SetMaxSize((max(dx,700),850))
    FPdlg.SetSize((max(dx,px+20),min(750,dy+py+30))) # 20 for scroll bar, 30 for a bit of room at bottom

def XferFPAsettings(InpParms):
    '''convert Topas-type parameters to SI units for NIST and place in a dict sorted
    according to use in each convoluter

    :param dict InpParms: a dict with Topas-like parameters, as set in
      :func:`MakeTopasFPASizer`
    :returns: a nested dict with global parameters and those for each convolution
    '''
    # cleanup old stuff
    for key in "tube_tails","absorption","si_psd","displacement","receiver_slit":
        if key in NISTparms:
            del NISTparms[key]

    keys = list(InpParms['wave'].keys())
    source_wavelengths_m = 1.e-10 * np.array([InpParms['wave'][i] for i in keys])
    la = [InpParms['int'][i] for i in keys]

    if IBmono:  # kludge: apply mono_slit_attenuation since it is not part of NIST FPA code
        norm = [InpParms['mono_slit_attenuation'] if
                    (1.5443 < InpParms['wave'][i] < 1.5447) else 1.
                for i in keys]
        source_intensities = norm * np.array(la)/max(la)
    else:
        source_intensities = np.array(la)/max(la)
    source_lor_widths_m = 1.e-10 * 1.e-3 * np.array([InpParms['lwidth'][i] for i in keys])
    source_gauss_widths_m = 1.e-10 * 1.e-3 * np.array([0.001 for i in keys])

    NISTparms["emission"] = {'emiss_wavelengths' : source_wavelengths_m,
                'emiss_intensities' : source_intensities,
                'emiss_gauss_widths' : source_gauss_widths_m,
                'emiss_lor_widths' : source_lor_widths_m,
                'crystallite_size_gauss' : 1.e-9 * InpParms.get('Size_G',1e6),
                'crystallite_size_lor' : 1.e-9 * InpParms.get('Size_L',1e6)}
    if IBmono:
        NISTparms["emission"]['a_mono'] = InpParms['src_mono_mm'] * 10**-3
        NISTparms["emission"]['b_mono'] = InpParms['focus_mono_mm'] * 10**-3
        NISTparms["emission"]['ibm_source_width'] = InpParms['mono_src_proj_mn'] * 10**-6
        for i in  ('passband_mistune','passband_shoulder','two_theta_mono'):
            NISTparms["emission"][i] = InpParms[i]
    elif InpParms.get('source_width', 0) > 0 and InpParms.get(
            'tube-tails_rel-I',0) > 0:
        NISTparms["tube_tails"] = {
            'main_width' : 1e-3 * InpParms.get('source_width', 0.),
            'tail_left' : -1e-3 * InpParms.get('tube-tails_L-tail',0.),
            'tail_right' : 1e-3 * InpParms.get('tube-tails_R-tail',0.),
            'tail_intens' : InpParms.get('tube-tails_rel-I',0.),}

    if InpParms['filament_length'] == InpParms['receiving_slit_length']: # workaround:
        InpParms['receiving_slit_length'] *= 1.00001 # avoid bug when slit lengths are identical
    NISTparms["axial"]  = {
            'axDiv':"full", 'slit_length_source' : 1e-3*InpParms['filament_length'],
            'slit_length_target' : 1e-3*InpParms['receiving_slit_length'],
            'length_sample' : 1e-3 * InpParms['sample_length'],
            'n_integral_points' : 10,
            'angI_deg' : InpParms['soller_angle'],
            'angD_deg': InpParms['soller_angle']
            }
    if InpParms.get('LAC_cm',0) > 0:
        NISTparms["absorption"] = {
            'absorption_coefficient': InpParms['LAC_cm']*100, #like LaB6, in m^(-1)
            'sample_thickness': 1e-3 * InpParms['sample_thickness'],
            }

    if InpParms.get('SiPSD_th2_angular_range',0) > 0 and DetMode == 'BBPSD':
        PSDdetector_length_mm=np.arcsin(np.pi*InpParms['SiPSD_th2_angular_range']/180.
                                            )*InpParms['Rs'] # mm
        NISTparms["si_psd"] = {
            'si_psd_window_bounds': (0.,PSDdetector_length_mm/1000.)
            }

    if InpParms.get('Specimen_Displacement'):
        NISTparms["displacement"] = {'specimen_displacement': 1e-3 * InpParms['Specimen_Displacement']}

    if InpParms.get('receiving_slit_width'):
        NISTparms["receiver_slit"] = {'slit_width':1e-3*InpParms['receiving_slit_width']}

    # set Global parameters
    max_wavelength = source_wavelengths_m[np.argmax(source_intensities)]
    NISTparms[""] = {
        'equatorial_divergence_deg' : InpParms['divergence'],
        'dominant_wavelength' : max_wavelength,
        'diffractometer_radius' : 1e-3* InpParms['Rs'],
        'oversampling' : InpParms['convolution_steps'],
        }

def setupFPAcalc():
    '''Create a peak profile object using the NIST XRD Fundamental
    Parameters Code.

    :returns: a profile object that can provide information on
      each convolution or compute the composite peak shape.
    '''
    if IBmono:
        p=FP.FP_windowed(anglemode="twotheta",
                    output_gaussian_smoother_bins_sigma=1.0,
                    oversampling=NISTparms.get('oversampling',10))
    else:
        p=FP.FP_profile(anglemode="twotheta",
                    output_gaussian_smoother_bins_sigma=1.0,
                    oversampling=NISTparms.get('oversampling',10))

    p.debug_cache=False
    #set parameters for each convolver
    for key in NISTparms:
        if key:
            p.set_parameters(convolver=key,**NISTparms[key])
        else:
            p.set_parameters(**NISTparms[key])
    return p

def doFPAcalc(NISTpk,ttArr,twotheta,calcwid,step):
    '''Compute a single peak using a NIST profile object

    :param object NISTpk: a peak profile computational object from the
      NIST XRD Fundamental Parameters Code, typically established from
      a call to :func:`SetupFPAcalc`
    :param np.Array ttArr: an evenly-spaced grid of two-theta points (degrees)
    :param float twotheta: nominal center of peak (degrees)
    :param float calcwid: width to perform convolution (degrees)
    :param float step: step size
    '''
    # find closest point to twotheta (may be outside limits of the array)
    center_bin_idx=min(ttArr.searchsorted(twotheta),len(ttArr)-1)
    NISTpk.set_optimized_window(twotheta_exact_bin_spacing_deg=step,
                twotheta_window_center_deg=ttArr[center_bin_idx],
                twotheta_approx_window_fullwidth_deg=calcwid,
                )
    NISTpk.set_parameters(twotheta0_deg=twotheta)
    return center_bin_idx,NISTpk.compute_line_profile()

def MakeSimSizer(G2frame, dlg):
    '''Create a GUI to get simulation with parameters for Fundamental
    Parameters fitting.

    :param wx.Window dlg: Frame or Dialog where GUI will appear

    :returns: a sizer with the GUI controls

    '''
    def _onOK(event):
        msg = ''
        if simParms['minTT']-simParms['calcwid']/1.5 < 0.1:
            msg += 'First peak minus half the calc width is too low'
        if simParms['maxTT']+simParms['calcwid']/1.5 > 175:
            if msg: msg += '\n'
            msg += 'Last peak plus half the calc width is too high'
        if simParms['npeaks'] < 8:
            if msg: msg += '\n'
            msg += 'At least 8 peaks are needed'
        if msg:
            G2G.G2MessageBox(dlg,msg,'Bad input, try again')
            return
        # compute "obs" pattern
        ttArr = np.arange(max(0.5,
                            simParms['minTT']-simParms['calcwid']/1.5),
                            simParms['maxTT']+simParms['calcwid']/1.5,
                            simParms['step'])
        intArr = np.zeros_like(ttArr)
        peaklist = np.linspace(simParms['minTT'],simParms['maxTT'],
                               simParms['npeaks'],endpoint=True)
        peakSpacing = (peaklist[-1]-peaklist[0])/(len(peaklist)-1)
        NISTpk = setupFPAcalc()
        minPtsHM = len(intArr)  # initialize points above half-max
        maxPtsHM = 0
        for num,twoth_peak in enumerate(peaklist):
            try:
                center_bin_idx,peakObj = doFPAcalc(
                    NISTpk,ttArr,twoth_peak,simParms['calcwid'],
                    simParms['step'])
            except Exception as err:
                if msg: msg += '\n'
                msg = "Error computing convolution, revise input. Error =\n"+str(err)
                continue
            if num == 0: G2plt.PlotFPAconvolutors(G2frame,NISTpk)
            pkMax = peakObj.peak.max()
            pkPts = len(peakObj.peak)
            minPtsHM = min(minPtsHM,sum(peakObj.peak >= 0.5*pkMax)) # points above half-max
            maxPtsHM = max(maxPtsHM,sum(peakObj.peak >= 0.5*pkMax)) # points above half-max
            startInd = center_bin_idx-(pkPts//2) #this should be the aligned start of the new data
            # scale peak so max I=10,000 and add into intensity array
            if startInd < 0:
                intArr[:startInd+pkPts] += 10000 * peakObj.peak[-startInd:]/pkMax
            elif startInd > len(intArr):
                break
            elif startInd+pkPts >= len(intArr):
                offset = pkPts - len( intArr[startInd:] )
                intArr[startInd:startInd+pkPts-offset] += 10000 * peakObj.peak[:-offset]/pkMax
            else:
                intArr[startInd:startInd+pkPts] += 10000 * peakObj.peak/pkMax
        # check if peaks are too closely spaced
        if maxPtsHM*simParms['step'] > peakSpacing/4:
            if msg: msg += '\n'
            msg += 'Maximum FWHM ({}) is too large compared to the peak spacing ({}). Decrease number of peaks or increase data range.'.format(
                maxPtsHM*simParms['step'], peakSpacing)
        # check if too few points across Hmax
        if minPtsHM < 10:
            if msg: msg += '\n'
            msg += 'There are only {} points above the half-max. 10 are needed. Dropping step size.'.format(minPtsHM)
            simParms['step'] *= 0.5
        if msg:
            G2G.G2MessageBox(dlg,msg,'Bad input, try again')
            wx.CallAfter(MakeSimSizer,G2frame, dlg)
            return
        # pattern has been computed successfully
        dlg.Destroy()
        wx.CallAfter(FitFPApeaks,ttArr, intArr, peaklist, maxPtsHM) # do peakfit outside event callback

    def FitFPApeaks(ttArr, intArr, peaklist, maxPtsHM):
        '''Perform a peak fit to the FP simulated pattern
        '''
        pgbar = wx.ProgressDialog('FPA Simulation','Starting FPA simulation',100,
            parent=G2frame,
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE #|wx.PD_CAN_ABORT
            )
        pgbar.Raise()
        wx.BeginBusyCursor()
        # pick out one or two most intense wavelengths
        ints = list(NISTparms['emission']['emiss_intensities'])
        Lam1 = NISTparms['emission']['emiss_wavelengths'][np.argmax(ints)]*1e10
        if 'two_theta_mono' in NISTparms['emission']: # is there an IBM?
            Lam2 = None   # Yes, ~monochromatic
        else:  # get lambda 2
            ints[np.argmax(ints)] = -1
            Lam2 = NISTparms['emission']['emiss_wavelengths'][np.argmax(ints)]*1e10
        histId = G2frame.AddSimulatedPowder(ttArr,intArr,
                                       'NIST Fundamental Parameters simulation',
                                       Lam1,Lam2)
        controls = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Controls'))
        controldat = controls.get('data',
                            {'deriv type':'analytic','min dM/M':0.001,})  #fil
        Parms,Parms2 = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,histId,'Instrument Parameters'))
        peakData = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,histId,'Peak List'))
        # set background to 0 with one term = 0; disable refinement
        bkg1,bkg2 = bkg = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,histId,'Background'))
        bkg1[1]=False
        bkg1[2]=0
        bkg1[3]=0.0
        limits = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,histId,'Limits'))
        # approximate asym correction
        try:
            Parms['SH/L'][1] = 0.25 * (
                NISTparms['axial']['length_sample']+
                NISTparms['axial']['slit_length_source']
                        ) / NISTparms['']['diffractometer_radius']
        except:
            pass
        pgbar.Update(5,newmsg='Creating peak list')
        pgbar.Raise()
        for pos in peaklist:
            i = ttArr.searchsorted(pos)
            area = sum(intArr[max(0,i-maxPtsHM):min(len(intArr),i+maxPtsHM)])
            peakData['peaks'].append(G2mth.setPeakparms(Parms,Parms2,pos,area))
        pgbar.Update(10,newmsg='Refining peak positions')
        histData = G2frame.GPXtree.GetItemPyData(histId)
        # refine peak positions only
        bxye = np.zeros(len(histData[1][1]))
        peakData['sigDict'] = G2pwd.DoPeakFit('LSQ',peakData['peaks'],
            bkg,limits[1],Parms,Parms2,histData[1],bxye,[],False,controldat)[0]
        pgbar.Update(20,newmsg='Refining peak positions && areas')
        # refine peak areas as well
        for pk in peakData['peaks']:
            pk[1] = True
        peakData['sigDict'] = G2pwd.DoPeakFit('LSQ',peakData['peaks'],
            bkg,limits[1],Parms,Parms2,histData[1],bxye,[],False,controldat)[0]
        pgbar.Update(40,newmsg='Refining profile function')
        # refine profile function
        for p in ('U', 'V', 'W', 'X', 'Y'):
            Parms[p][2] = True
        peakData['sigDict'] = G2pwd.DoPeakFit('LSQ',peakData['peaks'],
            bkg,limits[1],Parms,Parms2,histData[1],bxye,[],False,controldat)[0]
        pgbar.Update(70,newmsg='Refining profile function && asymmetry')
        # add in asymmetry
        Parms['SH/L'][2] = True
        peakData['sigDict'] = G2pwd.DoPeakFit('LSQ',peakData['peaks'],
            bkg,limits[1],Parms,Parms2,histData[1],bxye,[],False,controldat)[0]
        pgbar.Update(100,newmsg='Done')
        # reset "initial" profile
        for p in Parms:
            if len(Parms[p]) == 3:
                Parms[p][0] = Parms[p][1]
                Parms[p][2] = False
        pgbar.Destroy()
        wx.EndBusyCursor()
        # save Iparms
        pth = G2G.GetExportPath(G2frame)
        fldlg = wx.FileDialog(G2frame, 'Set name to save GSAS-II instrument parameters file', pth, '',
            'instrument parameter files (*.instprm)|*.instprm',wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if fldlg.ShowModal() == wx.ID_OK:
                filename = fldlg.GetPath()
                # make sure extension is .instprm
                filename = os.path.splitext(filename)[0]+'.instprm'
                File = open(filename,'w')
                G2fil.WriteInstprm(File, Parms)
                File.close()
                print ('Instrument parameters saved to: '+filename)
        finally:
            fldlg.Destroy()
        #GSASIIpath.IPyBreak()

    def _onClose(event):
        dlg.Destroy()
    def SetButtonStatus(done=False):
        OKbtn.Enable(bool(NISTparms))
        #saveBtn.Enable(bool(NISTparms))
        if done: _onOK(None)
    def _onSetFPA(event):
        # Create a non-modal dialog for Topas-style FP input.
        FPdlg = wx.Dialog(dlg,wx.ID_ANY,'FPA parameters',
                style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
        MakeTopasFPASizer(G2frame,FPdlg,SetButtonStatus)
        FPdlg.CenterOnParent()
        FPdlg.Raise()
        FPdlg.Show()
    def _onSaveFPA(event):
        filename = G2G.askSaveFile(G2frame,'','.NISTfpa',
                                       'dict of NIST FPA values',dlg)
        writeNIST(filename)

    def _onReadFPA(event):
        filename = G2G.GetImportFile(G2frame,
                message='Read file with dict of values for NIST Fundamental Parameters',
                parent=dlg,
                wildcard='dict of NIST FPA values|*.NISTfpa')
        if not filename: return
        if not filename[0]: return
        try:
            txt = open(filename[0],'r').read()
            NISTparms.clear()
            d = eval(txt)
            NISTparms.update(d)
        except Exception as err:
            G2G.G2MessageBox(dlg,
                    u'Error reading file {}:{}\n'.format(filename,err),
                    'Bad dict input')
        #GSASIIpath.IPyBreak()
        SetButtonStatus()

    if dlg.GetSizer(): dlg.GetSizer().Clear(True)
    MainSizer = wx.BoxSizer(wx.VERTICAL)
    topSizer = wx.BoxSizer(wx.HORIZONTAL)
    topSizer.Add(wx.StaticText(dlg,wx.ID_ANY,
            'Fit Profile Parameters to Peaks from Fundamental Parameters',
            style=wx.ALIGN_CENTER),0)
    topSizer.Add((5,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(dlg,helpIndex='FPA'))
    MainSizer.Add(topSizer,0,wx.EXPAND)
    G2G.HorizontalLine(MainSizer,dlg)
    MainSizer.Add((5,5),0)
    prmSizer = wx.FlexGridSizer(cols=2,hgap=3,vgap=5)
    text = wx.StaticText(dlg,wx.ID_ANY,'value',style=wx.ALIGN_CENTER)
    text.SetBackgroundColour(wx.WHITE)
    text.SetForegroundColour(wx.BLACK)
    prmSizer.Add(text,0,wx.EXPAND)
    text = wx.StaticText(dlg,wx.ID_ANY,'explanation',style=wx.ALIGN_CENTER)
    text.SetBackgroundColour(wx.WHITE)
    text.SetForegroundColour(wx.BLACK)
    prmSizer.Add(text,0,wx.EXPAND)
    for key,defVal,text in (
            ('minTT',3.,'Location of first peak in 2theta (deg)'),
            ('maxTT',123.,'Location of last peak in 2theta (deg)'),
            ('step',0.01,'Pattern step size (deg 2theta)'),
            ('npeaks',13,'Number of peaks'),
            ('calcwid',2.,'Range to compute each peak (deg 2theta)'),
            ):
        if key not in simParms: simParms[key] = defVal
        ctrl = G2G.ValidatedTxtCtrl(dlg,simParms,key,size=(70,-1))
        prmSizer.Add(ctrl,1,wx.ALL|WACV,1)
        txt = wx.StaticText(dlg,wx.ID_ANY,text,size=(300,-1))
        txt.Wrap(280)
        prmSizer.Add(txt)
    MainSizer.Add(prmSizer)

    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = wx.Button(dlg, wx.ID_ANY,'Input FP vals')
    btnsizer.Add(btn)
    btn.Bind(wx.EVT_BUTTON,_onSetFPA)
    #saveBtn = wx.Button(dlg, wx.ID_ANY,'Save FPA dict')
    #btnsizer.Add(saveBtn)
    #saveBtn.Bind(wx.EVT_BUTTON,_onSaveFPA)
    readBtn = wx.Button(dlg, wx.ID_ANY,'Read FPA dict')
    btnsizer.Add(readBtn)
    readBtn.Bind(wx.EVT_BUTTON,_onReadFPA)
    MainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER, 0)
    MainSizer.Add((-1,4),1,wx.EXPAND,1)
    txt = wx.StaticText(dlg,wx.ID_ANY,'If you use this, please cite: '+
                            G2G.GetCite('Fundamental parameter fitting'),
                            size=(350,-1))
    txt.Wrap(340)
    MainSizer.Add(txt,0,wx.ALIGN_CENTER)
    btnsizer = wx.BoxSizer(wx.HORIZONTAL)
    OKbtn = wx.Button(dlg, wx.ID_OK)
    OKbtn.SetDefault()
    btnsizer.Add(OKbtn)
    Cbtn = wx.Button(dlg, wx.ID_CLOSE,"Cancel")
    btnsizer.Add(Cbtn)
    MainSizer.Add(btnsizer, 0, wx.ALIGN_CENTER, 0)
    MainSizer.Add((-1,4),1,wx.EXPAND,1)
    # bindings for close of window
    OKbtn.Bind(wx.EVT_BUTTON,_onOK)
    Cbtn.Bind(wx.EVT_BUTTON,_onClose)
    SetButtonStatus()
    dlg.SetSizer(MainSizer)
    MainSizer.Layout()
    MainSizer.Fit(dlg)
    dlg.SetMinSize(dlg.GetSize())
    dlg.SendSizeEvent()
    dlg.Raise()

def GetFPAInput(G2frame):
    dlg = wx.Dialog(G2frame,wx.ID_ANY,'FPA input',
            style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
    MakeSimSizer(G2frame,dlg)
    dlg.CenterOnParent()
    dlg.Show()
    return

if __name__ == "__main__":
    app = wx.PySimpleApp()
    GSASIIpath.InvokeDebugOpts()
    frm = wx.Frame(None) # create a frame
    frm.Show(True)
    frm.TutorialImportDir = '/tmp'
    size = wx.Size(700,600)
    frm.plotFrame = wx.Frame(None,-1,'GSASII Plots',size=size,
                    style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX)
    frm.G2plotNB = G2plt.G2PlotNoteBook(frm.plotFrame,G2frame=frm)
    frm.plotFrame.Show()

    GetFPAInput(frm)

    app.MainLoop()
