# -*- coding: utf-8 -*-
'''
'''
import os
import shutil
import struct as st
import numpy as np
try:
    import xmltodict as xml
except:
    xml = None
import GSASIIobj as G2obj
import GSASIIfiles as G2fil
import GSASIIpath
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
            if GSASIIpath.condaTest():
                msg += ' To fix this press "Install packages" button below'
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
        def XtractXML():
            '''Read the XML info into the GSAS-II reader structure.
            This code from Bob with no significant mods.

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
            Tsum8 = 0
            isum8 = 0
            while i < nSteps :
                if datano:
                    entry = data['RawData']['DataRoutes']['DataRoute'][datano]['Datum'][i].split(',')
                else:
                    entry = data['RawData']['DataRoutes']['DataRoute']['Datum'][i].split(',')
                x[i] = float(entry[2])
                y[i] = float(entry[4])*float(entry[0])/effTime
                i = i + 1
                # both entry 6 & 8 seem to have a temperature in the multi-scan file
                # these columns are not present in the single-scan file. Using 8 (at random)
                try:
                    Tsum8 += float(entry[6])
                    isum8 += 1
                except:
                    pass
            if isum8 > 0:
                self.Sample['Temperature'] = 273. + Tsum8/isum8
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
                readOK = XtractXML()
                if readOK:
                    break
                #print(f"Rejected {fil} (no points)")
            except Exception as msg:
                print(f"Read of {fil} failed:\n{msg}")
        else:
            self.errors = "Error: no initial file found"
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
