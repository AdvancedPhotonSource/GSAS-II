# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2017-04-12 15:12:45 -0500 (Wed, 12 Apr 2017) $
# $Author: vondreele $
# $Revision: 2777 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASIIscriptable.py $
# $Id: GSASIIIO.py 2777 2017-04-12 20:12:45Z vondreele $
########### SVN repository information ###################
"""

"""
import os.path as ospath
import sys
import cPickle
def LoadDictFromProjFile(ProjFile):
    ''''Read a GSAS-II project file and load items to dictionary
    :param: str ProjFile: GSAS-II project (name.gpx) full file name
    :returns: dict Project: representation of gpx file following the GSAS-II tree struture
        for each item: key = tree name (e.g. 'Controls','Restraints',etc.), data is dict
        data dict = {'data':item data whch may be list, dict or None,'subitems':subdata (if any)}
    :returns: list nameList: names of main tree entries & subentries used to reconstruct project file
    Example for fap.gpx:
    Project = {                 #NB:dict order is not tree order
        u'Phases':{'data':None,'fap':{phase dict}}, 
        u'PWDR FAP.XRA Bank 1':{'data':[histogram data list],'Comments':comments,'Limits':limits, etc}, 
        u'Rigid bodies':{'data': {rigid body dict}},
        u'Covariance':{'data':{covariance data dict}},
        u'Controls':{'data':{controls data dict}}, 
        u'Notebook':{'data':[notebook list]}, 
        u'Restraints':{'data':{restraint data dict}}, 
        u'Constraints':{'data':{constraint data dict}}]}
    nameList = [                #NB: reproduces tree order
        [u'Notebook',],
        [u'Controls',], 
        [u'Covariance',], 
        [u'Constraints',],
        [u'Restraints',], 
        [u'Rigid bodies',], 
        [u'PWDR FAP.XRA Bank 1',
             u'Comments', 
             u'Limits', 
             u'Background', 
             u'Instrument Parameters',
             u'Sample Parameters', 
             u'Peak List', 
             u'Index Peak List', 
             u'Unit Cells List', 
             u'Reflection Lists'], 
        [u'Phases', 
             u'fap']]
    '''
    if not ospath.exists(ProjFile):
        print ('\n*** Error attempt to open project file that does not exist:\n   '+
            str(ProjFile))
        return
    file = open(ProjFile,'rb')
    print 'loading from file: ',ProjFile
    Project = {}
    nameList = []
    try:
        while True:
            try:
                data = cPickle.load(file)
            except EOFError:
                break
            datum = data[0]
            Project[datum[0]] = {'data':datum[1]}
            nameList.append([datum[0],])
            for datus in data[1:]:
                Project[datum[0]][datus[0]] = datus[1]
                nameList[-1].append(datus[0])
        file.close()
        print('project load successful')
    except:
        print("Error reading file "+str(ProjFile)+". This is not a GSAS-II .gpx file")
        return None
    return Project,nameList
    
def SaveDictToProjFile(Project,nameList,ProjFile):
    'Save a GSAS-II project file from dictionary/nameList created by LoadDictFromProjFile'
    file = open(ProjFile,'wb')
    print 'save to file: ',ProjFile
    for name in nameList:
        print name
        data = []
        item = Project[name[0]]
        print item
        data.append([name[0],item['data']])
        for item2 in name[1:]:
            data.append([item2,item[item2]])
        cPickle.dump(data,file,1)
    file.close()
    print('project save successful')

def main():
    'Needs a doc string'
    arg = sys.argv
    print arg
    if len(arg) > 1:
        GPXfile = arg[1]
        if not ospath.exists(GPXfile):
            print 'ERROR - ',GPXfile," doesn't exist!"
            exit()
        Project,nameList = LoadDictFromProjFile(GPXfile)
        SaveDictToProjFile(Project,nameList,'testout.gpx')
    else:
        print 'ERROR - missing filename'
        exit()
    print "Done"
         
if __name__ == '__main__':
    main()
