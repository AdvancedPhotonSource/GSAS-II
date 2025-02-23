# -*- coding: utf-8 -*-
'''
'''
import os
import shutil
import numpy as np
from .. import GSASIIpath
try:
    import xmltodict as xml
except Exception as msg:
    #if GSASIIpath.GetConfigValue('debug'): print(f'Debug: xmltodict error = {msg}')
    xml = None
from .. import GSASIIobj as G2obj

class brml_ReaderClass(G2obj.ImportPowderData):
    'Routines to import powder data from a zip Bruker .brml file'
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to self-reference
            extensionlist=('.brml',),
            strictExtension=True,
            formatName = 'Bruker brml',
            longFormatName = 'Bruker .brml powder data file'
            )
        if xml is None:
            self.UseReader = False
            msg = 'Bruker .brml Reader skipped because xmltodict module is not installed.'
            from .. import GSASIIfiles as G2fil
            G2fil.ImportErrorMsg(msg,{'Bruker .brml Importer':['xmltodict']})
        self.scriptable = True

    def ContentsValidator(self, filename):
        '''Validate by testing if the file can be opened by zip
        and if so, does it contain at least one file named "RawData*.xml"
        '''
        if xml is None:
            return False
        try:
            import zipfile as ZF
            with ZF.ZipFile(filename, 'r') as zipObj:
                for fil in zipObj.namelist():
                    if 'RawData' in fil and '.xml' in fil:
                        return True
                return False
        except:
            return False

    def Reader(self,filename, ParentFrame=None, **kwarg):
        'Read a Bruker brml file'
        def XtractXMLScan():
            '''Read the XML info into the GSAS-II reader structure.
            This structure seems to be for where the detector is scanned.
            Code from Bob with some modifications to read a wider
            assortment of files and sample temperature

            :returns: True if read suceeds. May also throw an exception
              on failure
            '''
            self.idstring = f'{os.path.basename(filename)} {os.path.basename(fil)}'
            self.powderentry[0] = filename
            self.powderentry[2] = filNum
            self.comments = []
            datano = 1
            try:
                nSteps = int(data['RawData']['DataRoutes']['DataRoute'][datano]['ScanInformation']['MeasurementPoints'])
            except KeyError:
                datano = 0
                nSteps = int(data['RawData']['DataRoutes']['DataRoute']['ScanInformation']['MeasurementPoints'])
            if nSteps <= 10: return False  # too short
            x = np.zeros(nSteps, dtype=float)
            y = np.zeros(nSteps, dtype=float)
            w = np.zeros(nSteps, dtype=float)

            if datano:
                effTime = float(data['RawData']['DataRoutes']['DataRoute'][datano]['ScanInformation']['TimePerStepEffective'])
            else:
                effTime = float(data['RawData']['DataRoutes']['DataRoute']['ScanInformation']['TimePerStepEffective'])

            # Extract 2-theta angle and counts from the XML document
            i=0
            T = {}
            for j in (5,6,8): # columns with temperature?
                T[f'{j}sum'] = 0
                T[f'{j}c'] = 0
            # data appears to be in column 4 most of the time
            # but at least in one file, it is in column 3
            y3 = []
            while i < nSteps :
                if datano:
                    entry = data['RawData']['DataRoutes']['DataRoute'][datano]['Datum'][i].split(',')
                else:
                    entry = data['RawData']['DataRoutes']['DataRoute']['Datum'][i].split(',')
                x[i] = float(entry[2])
                y[i] = float(entry[4])*float(entry[0])/effTime
                y3.append(float(entry[3]))
                i = i + 1
                # both entry 6 & 8 seem to have a temperature in the multi-scan file
                # these columns are not present in the single-scan file.
                # Try column 8 first
                for j in (5,6,8): # columns with temperature?
                    try:
                        t = float(entry[j])
                        if t > 0:
                            T[f'{j}sum'] += t
                            T[f'{j}c'] += 1
                    except:
                        pass
            #breakpoint()
            try:  # is there some range in col 4 values?
                if abs(1. - min(y)/max(y)) < 1e-4: raise Exception
            except:
                y = np.array(y3)
            w = np.where(y>0,1/y,0.)
            for j in (8,6,5): # columns with temperature?
                if T[f'{j}c'] > 0:
                    self.Sample['Temperature'] = 273. + T[f'{j}sum']/T[f'{j}c']
                    break
            self.powderdata = [x,y,w,np.zeros(nSteps),np.zeros(nSteps),np.zeros(nSteps)]
            return True

        def XtractXMLNoscan():
            '''Read the XML info into the GSAS-II reader structure.
            This structure seems to be for where the detector is stationary.

            :returns: True if read suceeds. May also throw an exception
              on failure
            '''
            self.idstring = f'{os.path.basename(filename)} {os.path.basename(fil)}'
            self.powderentry[0] = filename
            self.powderentry[2] = filNum
            self.comments = []
            try:
                scandata = data['RawData']['DataRoutes']['DataRoute']['ScanInformation']['ScaleAxes']
                if not scandata: return
                if 'ScaleAxisInfo' not in scandata: return
            except:
                return False
            try:
                start = float(scandata['ScaleAxisInfo']['Start'])
                stop = float(scandata['ScaleAxisInfo']['Stop'])
                incr = float(scandata['ScaleAxisInfo']['Increment'])
                nSteps = int(0.5 + (stop-start)/incr)
                if nSteps <= 10: return False  # too short
                xyT = data['RawData']['DataRoutes']['DataRoute']['Datum'].split(',')
            except:
                return False
            try:
                self.Sample['Temperature'] = 273. + float(xyT[2])  # seems to be temperature
            except:
                return
            y = np.array([float(y) for y in xyT[3:3+nSteps]])
            x = np.array([start + i * incr for i in range(len(y))])
            w = np.where(y>0,1/y,0.)
            self.powderdata = [x,y,w,np.zeros(nSteps),np.zeros(nSteps),np.zeros(nSteps)]
            return True

        #### beginning of Reader
        if xml is None:
            return False
        self.buffer = kwarg.get('buffer',{})
        if 'filesFound' not in self.buffer or 'maxNum' not in self.buffer:
            # do a scan of zip to find RawData
            self.buffer['filesFound'] = []
            self.buffer['maxNum'] = -1
            try:
                import zipfile as ZF
                with ZF.ZipFile(filename, 'r') as zipObj:
                    for fil in zipObj.namelist():
                        if 'RawData' in fil and '.xml' in fil:
                            self.buffer['filesFound'].append(fil)
                        try:
                            self.buffer['maxNum'] = max(
                                self.buffer['maxNum'],
                                int(fil.split('RawData')[1].split('.')[0]))
                        except:
                            pass
            except Exception as msg:
                self.errors = "Error with dir unzip:\n"+str(msg)
                return False
        filNum = kwarg.get('blocknum',1)-1
        # loop over the files in the zip structure until we find one we can read
        while filNum <= self.buffer['maxNum']:
            # find file(s) matching the number
            files = [i for i,j in enumerate(self.buffer['filesFound'])
                         if f'RawData{filNum}.xml' in j]
            if not files:
                filNum += 1
                continue
            indx = files[0]
            fil = self.buffer['filesFound'].pop(indx) # take the 1st located file
            try:
                import zipfile as ZF
                with ZF.ZipFile(filename, 'r') as zipObj:
                    zipObj.extract(fil)
                with open(fil) as fd:
                    data = dict(xml.parse(fd.read()))
            except Exception as msg:
                print(f"Error with unzip/xml extraction on {fil}:\n{msg}")
                continue
            finally:
                if os.path.exists(os.path.dirname(fil)):
                    shutil.rmtree(os.path.dirname(fil))
            try:
                readOK = XtractXMLNoscan()
                if readOK: break
                readOK = XtractXMLScan()
                if readOK: break
                #print(f"Rejected {fil} (no points)")
            except Exception as msg:
                print(f"Read of {fil} failed:\n{msg}")
        else:
            self.errors = "Error: no datasets found to read"
            return False
        try:
            import zipfile as ZF
            with ZF.ZipFile(filename, 'r') as zipObj:
                zipObj.extract(fil)
            with open(fil) as fd:
                data = dict(xml.parse(fd.read()))
            os.remove(fil)
            os.rmdir('Experiment0')
        except Exception as msg:
            self.errors = f"Error with unzip/read:\n{msg}"
            return False
        finally:
            if os.path.exists(os.path.dirname(fil)):
                shutil.rmtree(os.path.dirname(fil))
        # add some comments
        for key in 'TimeStampStarted','TimeStampFinished':
            self.comments.append(f'{key}: {data["RawData"].get(key,"?")}\n')
        for key in data['RawData'].get('Identifier',{}):
            val = data['RawData']['Identifier'][key]
            if type(val) is str:
                self.comments.append(f'{key}: {val}\n')
        # file read was successful, are there more files to read?
        self.repeatcount = filNum
        while filNum <= self.buffer['maxNum']:
            # find file(s) matching the number
            files = [i for i,j in enumerate(self.buffer['filesFound'])
                         if f'RawData{filNum}.xml' in j]
            if not files:
                filNum += 1
                continue
            self.repeat = True
            break
        else:
            self.repeat = False
        return True
