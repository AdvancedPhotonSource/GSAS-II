# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: $
# $Author: von dreele $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################

from __future__ import division, print_function
import os.path as ospath
import xml.etree.ElementTree as ET
import numpy as np
import GSASIIobj as G2obj
import GSASIIpath
sind = lambda x: np.sin(x*np.pi/180.)
GSASIIpath.SetVersionNumber("$Revision: $")
class Panalytical_ReaderClass(G2obj.ImportReflectometryData):
    '''Routines to import reflectivity data from a Pananalytical.xrdm (xml) file. 
    
    '''
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.xrdml','.xml'),
            strictExtension=True,
            formatName = 'Panalytical xrdml (xml)',
            longFormatName = 'Panalytical reflectivity data as *.xrdml'
            )
        self.scriptable = True
        self.vals = None
        self.stepsize = None
        self.skip = 0
        self.root = None

    # Validate the contents -- make sure we only have valid lines and set
    # values we will need for later read.
    def ContentsValidator(self, filename):
        fp = open(filename,'r')
        self.vals = None
        self.stepsize = None
        fp.seek(0)
        try:
            self.root = ET.parse(fp).getroot()
            tag = self.root.tag
            tag = tag.split('}')[0]+'}'
            self.root.find(tag+'comment')
        except:
            self.errors = 'Bad xml file'
            fp.close()
            return False
        fp.close()
        return True
            
    def Reader(self,filename, ParentFrame=None, **kwarg):
        'Read a Panalytical .xrdml (.xml) file; already in self.root'
        blockNum = kwarg.get('blocknum',0)
        self.idstring = ospath.basename(filename) + ' Scan '+str(blockNum)
        x = []
        y = []
        w = []
        tag = self.root.tag
        tag = tag.split('}')[0]+'}'
        sample = self.root.find(tag+'sample')
        self.idstring = ospath.basename(filename) + ' Scan '+str(blockNum)
        blks = self.root.findall(tag+'xrdMeasurement')
        scans = []
        for data in blks:
            scans += data.findall(tag+'scan')
        data = self.root.find(tag+'xrdMeasurement')        
        wave = data.find(tag+'usedWavelength')
        incident = data.find(tag+'incidentBeamPath')
        radius = float(incident.find(tag+'radius').text)
        tube = incident.find(tag+'xRayTube')
        if len(scans) > 1:
            self.repeat = True
        if blockNum-1 == len(scans):
            self.repeat = False
            return False
        scan = scans[blockNum-1]
        header = scan.find(tag+'header')
        dataPoints = scan.find(tag+'dataPoints')
        self.comments.append('Gonio. radius=%.2f'%(radius))
        self.Sample['Gonio. radius'] = radius
        if sample.find(tag+'id').text:
            self.comments.append('Sample name='+sample.find(tag+'id').text)
        try:
            self.comments.append('Date/TimeStart='+header.find(tag+'startTimeStamp').text)
            self.comments.append('Date/TimeEnd='+header.find(tag+'endTimeStamp').text)
            self.comments.append('xray tube='+tube.attrib['name'])
        except AttributeError:
            pass
        self.comments.append('Ka1=%s'%(wave.find(tag+'kAlpha1').text))
        self.comments.append('Ka2=%s'%(wave.find(tag+'kAlpha2').text))
        self.comments.append('Ka2/Ka1=%s'%(wave.find(tag+'ratioKAlpha2KAlpha1').text))
        self.comments.append('Kb=%s'%(wave.find(tag+'kBeta').text))
        self.comments.append('Voltage='+tube.find(tag+'tension').text)
        self.comments.append('Current='+tube.find(tag+'current').text)
        limits = dataPoints.find(tag+'positions')
        startPos = float(limits.find(tag+'startPosition').text)
        endPos= float(limits.find(tag+'endPosition').text)
        for seclbl in 'intensities','counts':
            sec = dataPoints.find(tag+seclbl)
            if sec is None: continue
            y = np.fromstring(sec.text,sep=' ')
            break
        else:
            print('Panalytical read error: Intensities could not be located')
            return False            
        self.instdict['wave'] = float(wave.find(tag+'kAlpha1').text)
        self.instdict['type'] = 'RXC'
        self.reflectometryentry[0] = filename
        self.reflectometryentry[2] = blockNum 
        N = y.shape[0]
        x = np.linspace(startPos,endPos,N)
        x = 4.*np.pi*sind(x/2.)/self.instdict['wave']
        w = np.where(y>0,1./y,1.)
        self.reflectometrydata = [
            np.array(x), # x-axis values
            np.array(y), # powder pattern intensities
            np.array(w), # 1/sig(intensity)^2 values (weights)
            np.zeros(N), # calc. intensities (zero)
            np.zeros(N), # calc. background (zero)
            np.zeros(N), # obs-calc profiles
            ]
        conditions = scan.find(tag+'nonAmbientPoints')
        if conditions is not None:
            kind = conditions.attrib['type']
            if kind == 'Temperature':
                Temperature = float(conditions.find(tag+'nonAmbientValues').text.split()[-1])
                self.Sample['Temperature'] = Temperature
        return True
