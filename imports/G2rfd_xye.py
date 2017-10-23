# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2017-03-04 08:34:05 -0600 (Sat, 04 Mar 2017) $
# $Author: vondreele $
# $Revision: 2738 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/G2sad_xye.py $
# $Id: G2sad_xye.py 2738 2017-03-04 14:34:05Z vondreele $
########### SVN repository information ###################
'''
*Module G2rfd_xye: read reflectometry data*
------------------------------------------------

Routines to read in reflectometry data from an .xye type file, with
two-theta or Q steps. 

'''

from __future__ import division, print_function
import os.path as ospath
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 2738 $")
npasind = lambda x: 180.*np.arcsin(x)/np.pi
npsind = lambda x: np.sin(np.pi*x/180.)
fourpi = 4.0*np.pi

class txt_XRayReaderClass(G2obj.ImportReflectometryData):
    'Routines to import X-ray q REFD data from a .xrfd or .xdat file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.xrfd','.xdat'),
            strictExtension=False,
            formatName = 'q (A-1) step X-ray QRE data',
            longFormatName = 'q (A-1) stepped X-ray text data file in Q,R,E order; E optional'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid q-step file'
        fp = open(filename,'r')
        Ndata = 0
        for i,S in enumerate(fp):
            if '#' in S[0]:
                continue
            vals = S.split()
            if len(vals) >= 2:
                try:
                    data = [float(val) for val in vals]
                    Ndata += 1
                except ValueError:
                    pass
        fp.close()
        if not Ndata:     
            self.errors = 'No 2 or more column numeric data found'
            return False
        return True # no errors encountered

    def Reader(self,filename, ParentFrame=None, **unused):
        print ('Read a q-step text file')
        x = []
        y = []
        w = []
        sq = []
        wave = 1.5428   #Cuka default
        Temperature = 300
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if len(S) == 1:     #skip blank line
                continue
            if '=' in S:
                self.comments.append(S[:-1])
                if 'wave' in S.split('=')[0].lower():
                    try:
                        wave = float(S.split('=')[1])
                    except:
                        pass
                continue
            if '#' in S[0]:
                continue
            vals = S.split()
            if len(vals) >= 2:
                try:
                    data = [float(val) for val in vals]
                    x.append(float(data[0]))
                    f = float(data[1])
                    if f <= 0.0:
                        del x[-1]
                        continue
                    elif len(vals) > 2:
                        y.append(float(data[1]))
                        w.append(1.0/float(data[2])**2)
                        if len(vals) == 4:
                            sq.append(float(data[3]))
                        else:
                            sq.append(0.)
                    else:
                        y.append(float(data[1]))
                        w.append(1.0/(0.02*float(data[1]))**2)
                        sq.append(0.)
                except ValueError:
                    msg = 'Error in line '+str(i+1)
                    print (msg)
                    continue
        fp.close()
        N = len(x)
        for S in self.comments:
            if 'Temp' in S.split('=')[0]:
                try:
                    Temperature = float(S.split('=')[1])
                except:
                    pass
        self.instdict['wave'] = wave
        self.instdict['type'] = 'RXC'
        x = np.array(x)
        self.reflectometrydata = [
            x, # x-axis values q
            np.array(y), # small angle pattern intensities
            np.array(w), # 1/sig(intensity)^2 values (weights)
            np.zeros(N), # calc. intensities (zero)
            np.zeros(N), # obs-calc profiles
            np.array(sq), # fix bkg
            ]
        self.reflectometryentry[0] = filename
        self.reflectometryentry[2] = 1 # xye file only has one bank
        self.idstring = ospath.basename(filename)
        # scan comments for temperature
        self.Sample['Temperature'] = Temperature
        return True

class txt_NeutronReaderClass(G2obj.ImportReflectometryData):
    'Routines to import neutron q REFD data from a .nrfd or .ndat file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.nrfd','.ndat'),
            strictExtension=False,
            formatName = 'q (A-1) step neutron QRE data',
            longFormatName = 'q (A-1) stepped neutron text data file in Q,R,E order; E optional'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid q-step file'
        Ndata = 0
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if '#' in S[0]:
                continue
            vals = S.split()
            if len(vals) >= 2:
                try:
                    data = [float(val) for val in vals]
                    Ndata += 1
                except ValueError:
                    pass
        fp.close()
        if not Ndata:     
            self.errors = 'No 2 or more column numeric data found'
            return False
        return True # no errors encountered

    def Reader(self,filename,filepointer, ParentFrame=None, **unused):
        print ('Read a q-step text file')
        x = []
        y = []
        w = []
        sq = []
        wave = 1.5428   #Cuka default
        Temperature = 300
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if len(S) == 1:     #skip blank line
                continue
            if '=' in S:
                self.comments.append(S[:-1])
                if 'wave' in S.split('=')[0].lower():
                    try:
                        wave = float(S.split('=')[1])
                    except:
                        pass
                continue
            if '#' in S[0]:
                continue
            vals = S.split()
            if len(vals) >= 2:
                try:
                    data = [float(val) for val in vals]
                    x.append(float(data[0]))
                    f = float(data[1])
                    if f <= 0.0:
                        del x[-1]
                        continue
                    elif len(vals) > 2:
                        y.append(float(data[1]))
                        w.append(1.0/float(data[2])**2)
                        if len(vals) == 4:
                            sq.append(float(data[3]))
                        else:
                            sq.append(0.)
                    else:
                        y.append(float(data[1]))
                        w.append(1.0/(0.02*float(data[1]))**2)
                        sq.append(0.)
                except ValueError:
                    msg = 'Error in line '+str(i+1)
                    print (msg)
                    continue
        fp.close()
        N = len(x)
        for S in self.comments:
            if 'Temp' in S.split('=')[0]:
                try:
                    Temperature = float(S.split('=')[1])
                except:
                    pass
        self.instdict['wave'] = wave
        self.instdict['type'] = 'RNC'
        x = np.array(x)
        self.reflectometrydata = [
            x, # x-axis values q
            np.array(y), # small angle pattern intensities
            np.array(w), # 1/sig(intensity)^2 values (weights)
            np.zeros(N), # calc. intensities (zero)
            np.zeros(N), # obs-calc profiles
            np.array(sq), # Q FWHM
            ]
        self.reflectometryentry[0] = filename
        self.reflectometryentry[2] = 1 # xye file only has one bank
        self.idstring = ospath.basename(filename)
        # scan comments for temperature
        self.Sample['Temperature'] = Temperature
        return True

class txt_XRayThetaReaderClass(G2obj.ImportReflectometryData):
    'Routines to import X-ray theta REFD data from a .xtrfd or .xtdat file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.xtrfd','.xtdat'),
            strictExtension=False,
            formatName = 'theta step X-ray QRE data',
            longFormatName = 'theta stepped X-ray text data file in Q,R,E order; E optional'
            )

    # Validate the contents -- make sure we only have valid lines
    def ContentsValidator(self, filename):
        'Look through the file for expected types of lines in a valid q-step file'
        Ndata = 0
        self.wavelength = 0.
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if '#' in S[0]:
                if 'wavelength' in S[:-1].lower():
                    self.wavelength = float(S[:-1].split('=')[1])
                elif 'energy' in S[:-1].lower():
                    self.wavelength = 12.39842*1000./float(S[:-1].split('=')[1])
                continue
            vals = S.split()
            if len(vals) >= 2:
                try:
                    data = [float(val) for val in vals]
                    Ndata += 1
                except ValueError:
                    pass
        fp.close()
        if not Ndata:     
            self.errors = 'No 2 or more column numeric data found'
            return False
        elif not self.wavelength:
            self.errors = 'Missing wavelength or energy in header'
            return False
        return True # no errors encountered

    def Reader(self,filename, ParentFrame=None, **unused):
        print ('Read a q-step text file')
        x = []
        y = []
        w = []
        sq = []
        wave = self.wavelength
        Temperature = 300
        fp = open(filename,'r')
        for i,S in enumerate(fp):
            if len(S) == 1:     #skip blank line
                continue
            if '=' in S:
                self.comments.append(S[:-1])
                if 'wave' in S.split('=')[0].lower():
                    try:
                        wave = float(S.split('=')[1])
                    except:
                        pass
                continue
            if '#' in S[0]:
                continue
            vals = S.split()
            if len(vals) >= 2:
                try:
                    data = [float(val) for val in vals]
                    x.append(fourpi*npsind(float(data[0]))/wave)
                    f = float(data[1])
                    if f <= 0.0:
                        del x[-1]
                        continue
                    elif len(vals) > 2:
                        y.append(float(data[1]))
                        w.append(1.0/float(data[2])**2)
                        if len(vals) == 4:
                            sq.append(float(data[3]))
                        else:
                            sq.append(0.)
                    else:
                        y.append(float(data[1]))
                        w.append(1.0/(0.02*float(data[1]))**2)
                        sq.append(0.)
                except ValueError:
                    msg = 'Error in line '+str(i+1)
                    print (msg)
                    continue
        fp.close()
        N = len(x)
        for S in self.comments:
            if 'Temp' in S.split('=')[0]:
                try:
                    Temperature = float(S.split('=')[1])
                except:
                    pass
        self.instdict['wave'] = wave
        self.instdict['type'] = 'RXC'
        x = np.array(x)
        self.reflectometrydata = [
            x, # x-axis values q
            np.array(y), # small angle pattern intensities
            np.array(w), # 1/sig(intensity)^2 values (weights)
            np.zeros(N), # calc. intensities (zero)
            np.zeros(N), # obs-calc profiles
            np.array(sq), # fix bkg
            ]
        self.reflectometryentry[0] = filename
        self.reflectometryentry[2] = 1 # xye file only has one bank
        self.idstring = ospath.basename(filename)
        # scan comments for temperature
        self.Sample['Temperature'] = Temperature

        return True

