# -*- coding: utf-8 -*-
#GSASIIrestr - restraint GUI routines
'''Restraint GUI routines follow.
'''
from __future__ import division, print_function
import wx
import wx.grid as wg
import numpy as np
import numpy.ma as ma
import os.path
from . import GSASIIpath
from . import GSASIImath as G2mth
from . import GSASIIlattice as G2lat
from . import GSASIIspc as G2spc
from . import GSASIIdataGUI as G2gd
from . import GSASIIplot as G2plt
from . import GSASIIdata as G2data
from . import GSASIIctrlGUI as G2G
from . import GSASIIphsGUI as G2phsGUI
from . import GSASIIobj as G2obj
from . import GSASIIconstrGUI as G2cnstG
from . import GSASIIexprGUI as G2exG

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

try:
    VERY_LIGHT_GREY = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE)
    VERY_YELLOW = wx.Colour(255,255,0)
    WACV = wx.ALIGN_CENTER_VERTICAL
except:
    pass
TabSelectionIdDict = {}

################################################################################
#####  Restraints
################################################################################           
def GetSelectedRows(widget,lbl='edit',G2frame=None):
    '''Returns a list of selected rows. Rows can be selected, blocks of cells
    or individual cells can be selected. The column for selected cells is ignored.
    '''
    try:
        rows = widget.GetSelectedRows()
        if rows: return rows
            
        top = widget.GetSelectionBlockTopLeft()
        bot = widget.GetSelectionBlockBottomRight()
        if top and bot:
            rows = range(top[0][0],bot[0][0]+1)
            if rows: return rows

        rows = sorted(list(set([cell[0] for cell in widget.GetSelectedCells()])))
        if rows: return rows

        choices = ["{}: {}".format(widget.GetRowLabelValue(i),widget.GetCellValue(i,0)) 
                       for i in range(widget.GetNumberRows())]
        try:
            dlg = G2G.G2MultiChoiceDialog(G2frame,'Restraints to '+lbl,
                                                  'Select restraints',choices)
            if dlg.ShowModal() != wx.ID_OK: return
            return dlg.GetSelections()
        finally:
            dlg.Destroy()
    except:
        return
    
def UpdateRestraints(G2frame,data,phaseName):
    '''Respond to selection of the Restraints item on the
    data tree
    '''
#    global Pages
    
    def getMacroFile(macName):
        defDir = os.path.join(GSASIIpath.path2GSAS2,'inputs','GSASIImacros')
        if not os.path.exists(defDir): # patch 3/2024 for svn dir organization
            defDir = os.path.join(GSASIIpath.path2GSAS2,'GSASIImacros')
        if not os.path.exists(defDir):
            print('Warning: GSASIImacros directory not found')
            return
        dlg = wx.FileDialog(G2frame,message='Choose '+macName+' restraint macro file',
            defaultDir=defDir,defaultFile="",wildcard="GSAS-II macro file (*.mac)|*.mac",
            style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            macro = ''
            if dlg.ShowModal() == wx.ID_OK:
                macfile = dlg.GetPath()
                macro = open(macfile,'r')
                head = macro.readline()
                if macName not in head:
                    print (head)
                    print ('**** ERROR - wrong restraint macro file selected, try again ****')
                    macro = []
        finally:
            dlg.Destroy()
        return macro        #advanced past 1st line
        
    def getMOGULFile():
        colNums = [0,2,3,5,6,7] # location for these fields:
        # Type, Fragment, No. of hits, Query value, Mean, Std. dev.
        dlg = wx.FileDialog(G2frame,message='Choose MOGUL csv file',
            defaultDir='.',defaultFile="",wildcard="MOGUL csv file (*.csv)|*.csv",
            style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        try:
            mogul = ''
            if dlg.ShowModal() == wx.ID_OK:
                csvfile = dlg.GetPath()
                mogul = open(csvfile,'r')
                head = mogul.readline()
                if 'Type' not in head:
                    print ('Note: header line is\n',head)
                    print ('**** ERROR - file selected is not a MOGUL csv file, try again ****')
                    mogul = []
                else:
                    for i,k in enumerate(('Type','Fragment',
                            'No. of hits','Query value','Mean','Std. dev.')):
                        try:
                            colNums[i] = head.split(',').index(k)
                        except ValueError:
                            pass
        finally:
            dlg.Destroy()
        return mogul,colNums        #file pointer advanced past 1st line & col pointers
        
    def OnPlotAARestraint(event):
        page = G2frame.restrBook.GetSelection()
        if 'Torsion' in G2frame.restrBook.GetPageText(page):
            torNames = []
            torNames += list(restrData['Torsion']['Coeff'].keys())
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Torsion data',torNames)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    torName = torNames[dlg.GetSelection()]
                    torsion = G2data.torsionDist[torName]
                    torCoeff = restrData['Torsion']['Coeff'][torName]
                    torList = restrData['Torsion']['Torsions']
                    Names = []
                    Angles = []
                    for i,[indx,ops,cofName,esd] in enumerate(torList):
                        if cofName == torName:
                            atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                            name = '('+atoms[2][1]+atoms[2][0].strip()+atoms[2][2]+')'
                            for atom in atoms:
                                name += '  '+atom[3]
                            XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                            angle = G2mth.getRestTorsion(XYZ,Amat)
                            Angles.append(angle)
                            Names.append(name) 
                    G2plt.PlotTorsion(G2frame,phaseName,torsion,torName,Names,np.array(Angles),torCoeff)
            finally:
                dlg.Destroy()
            
        elif 'Rama' in G2frame.restrBook.GetPageText(page):
            ramaNames = ['All',]
            ramaNames += list(restrData['Rama']['Coeff'].keys())
            dlg = wx.SingleChoiceDialog(G2frame,'Select','Ramachandran data',ramaNames)
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    ramaName = ramaNames[dlg.GetSelection()]
                    rama = G2data.ramachandranDist[ramaName]
                    ramaCoeff = []
                    if ramaName != 'All':
                        ramaCoeff = restrData['Rama']['Coeff'][ramaName]
                    ramaList = restrData['Rama']['Ramas']
                    Names = []
                    PhiPsi = []
                    for i,[indx,ops,cofName,esd] in enumerate(ramaList):
                        if cofName == ramaName or (ramaName == 'All' and '-1' in cofName):
                            atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                            name = '('+atoms[3][1]+atoms[3][0].strip()+atoms[3][2]+')'
                            for atom in atoms:
                                name += '  '+atom[3]
                            XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                            phi,psi = G2mth.getRestRama(XYZ,Amat)
                            PhiPsi.append([phi,psi])
                            Names.append(name) 
                    G2plt.PlotRama(G2frame,phaseName,rama,ramaName,Names,np.array(PhiPsi),ramaCoeff)
            finally:
                dlg.Destroy()

    def SetupParmDict(G2frame):
        '''Creates a parameter dict with variable names as keys and 
        numerical values (only)
        '''
        G2cnstG.CheckAllScalePhaseFractions(G2frame,refine=False)
        try:
            parmDict,varyList = G2frame.MakeLSParmDict()
        except:
            print('Error retrieving parameters')
            return {}
        return {i:parmDict[i][0] for i in parmDict}
        
    def OnAddRestraint(event):
        '''Adds a restraint depending on which tab is currently displayed'''
        page = G2frame.restrBook.GetSelection()
        if 'Bond' in G2frame.restrBook.GetPageText(page):
            AddBondRestraint(restrData['Bond'])
        elif 'Angle' in G2frame.restrBook.GetPageText(page):
            AddAngleRestraint(restrData['Angle'])
        elif 'Plane' in G2frame.restrBook.GetPageText(page):
            AddPlaneRestraint(restrData['Plane'])
        elif 'Chem' in G2frame.restrBook.GetPageText(page):
            AddChemCompRestraint(restrData['ChemComp'])
        elif 'Texture' in G2frame.restrBook.GetPageText(page):
            AddTextureRestraint(restrData['Texture'])
        elif 'Moments' in G2frame.restrBook.GetPageText(page):
            AddMomentRestraint(restrData['Moments'])
        elif 'General' in G2frame.restrBook.GetPageText(page):
            parmDict = SetupParmDict(G2frame)
            dlg = G2exG.ExpressionDialog(G2frame,parmDict,
                           header="Create a restraint expression",
                           fit=False,wildCard=G2frame.testSeqRefineMode())
            restobj = dlg.Show(True)
            if restobj:
                restrData['General']['General'].append([restobj,0.0,1.0])
                wx.CallAfter(UpdateGeneralRestr,restrData['General'])
            
    def OnAddAARestraint(event):
        page = G2frame.restrBook.GetSelection()
        if 'Bond' in G2frame.restrBook.GetPageText(page):
            AddAABondRestraint(restrData['Bond'])
        elif 'Angle' in G2frame.restrBook.GetPageText(page):
            AddAAAngleRestraint(restrData['Angle'])
        elif 'Plane' in G2frame.restrBook.GetPageText(page):
            AddAAPlaneRestraint(restrData['Plane'])
        elif 'Chiral' in G2frame.restrBook.GetPageText(page):
            AddAAChiralRestraint(restrData['Chiral'])
        elif 'Torsion' in G2frame.restrBook.GetPageText(page):
            AddAATorsionRestraint(restrData['Torsion'])
        elif 'Rama' in G2frame.restrBook.GetPageText(page):
            AddAARamaRestraint(restrData['Rama'])
            
    def OnUseMogul(event):
        page = G2frame.restrBook.GetSelection()
        if 'Bond' in G2frame.restrBook.GetPageText(page):
            AddMogulBondRestraint(restrData['Bond'])
        elif 'Angle' in G2frame.restrBook.GetPageText(page):
            AddMogulAngleRestraint(restrData['Angle'])
            
    def makeChains(Names,Ids):
        Chains = {}
        atoms = zip(Names,Ids)
        for name,Id in atoms:
            items = name.split(' ',2)
            rnum,res = items[0].split(':')
            rnum = int(rnum)
            if items[1] not in Chains:
                Residues = {}
                Chains[items[1]] = Residues
            if rnum not in Residues:
                Residues[rnum] = [[],[]]
            if items[2][3] in [' ','A']:
                Residues[rnum][0].append([res,items[2],Id])
            if items[2][3] in [' ','B']:
                Residues[rnum][1].append([res,items[2],Id])
        return Chains
        
    def AddBondRestraint(bondRestData):
        Lists = {'origin':[],'target':[]}
        for listName in ['origin','target']:
            dlg = G2G.G2MultiChoiceDialog(G2frame,'Bond restraint '+listName+' for '+General['Name'],
                    'Select bond restraint '+listName+' atoms',Names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                for x in sel:
                    if 'all' in Names[x]:
                        allType = Types[x]
                        for name,Type,coords,Id in zip(Names,Types,Coords,Ids):
                            if Type == allType and 'all' not in name:
                                Lists[listName].append([Id,Type,coords])
                    else:
                        Lists[listName].append([Ids[x],Types[x],Coords[x],])
            else:
                return
        if len(Lists['origin']) and len(Lists['target']):
            bond = 1.54
            dlg = G2G.SingleFloatDialog(G2frame,'Distance','Enter restraint distance for bond',bond,[0.01,4.],'%.4f')
            if dlg.ShowModal() == wx.ID_OK:
                bond = dlg.GetValue()
            dlg.Destroy()
        Factor = bondRestData['Range']
        dlg = wx.ProgressDialog("Generating bond restraints","Processed origin atoms",len(Lists['origin']), 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME)
        try:
            bondlst = G2mth.searchBondRestr(Lists['origin'],Lists['target'],
                bond,Factor,General['Type'],SGData,Amat,0.01,dlg)
            for newBond in bondlst:
                if newBond not in bondRestData['Bonds']:
                    bondRestData['Bonds'].append(newBond)
                
        finally:
            dlg.Destroy()
        UpdateBondRestr(bondRestData)                

    def AddAABondRestraint(bondRestData):
        macro = getMacroFile('bond')
        if not macro:
            return
        macStr = macro.readline()
        atoms = zip(Names,Coords,Ids)
        Factor = bondRestData['Range']
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Bond']['wtFactor'] = float(items[1])
            elif 'S' in items[0]:
                oIds = []
                oCoords = []
                oDis = []
                tIds = []
                tCoords = []
                tDis = []
                res = items[1]
                dist = float(items[2])
                esd = float(items[3])
                oAtm,tAtm = items[4:6]
                for Name,coords,Id in atoms:
                    names = Name.split(' ',2)
                    if res == '*' or res in names[0]:
                        if oAtm.ljust(3) == names[2][:3]:
                            oIds.append(Id)
                            oCoords.append(np.array(coords))
                            oDis.append(names[2][3])
                        if tAtm.ljust(3) == names[2][:3]:
                            tIds.append(Id)
                            tCoords.append(np.array(coords))
                            tDis.append(names[2][3])
                for i,[oId,oCoord,odis] in enumerate(zip(oIds,oCoords,oDis)):
                    for tId,tCoord,tdis in list(zip(tIds,tCoords,tDis))[i:]:
                        if odis+tdis in ['AB','BA']:
                            continue
                        obsd = np.sqrt(np.sum(np.inner(Amat,tCoord-oCoord)**2))
                        if dist/Factor < obsd < dist*Factor:
                            newBond = [[oId,tId],['1','1'],dist,esd]
                            if newBond not in bondRestData['Bonds']:
                                bondRestData['Bonds'].append(newBond)              
            macStr = macro.readline()
        macro.close()
        print(' Found %d bond restraints'%len(bondRestData['Bonds']))
        UpdateBondRestr(bondRestData)

    def AddMogulBondRestraint(bondRestData):
        mogul,colNums = getMOGULFile()
        badNames = []
        badCount = 0
        for line in mogul:
            items = line.split(',')
            if 'bond' == items[colNums[0]]:
                oName,tName = items[colNums[1]].split()
                try:
                    oInd = Names.index(oName)
                except:
                    badCount += 1
                    badNames.append(oName)
                    continue
                try:
                    tInd = Names.index(tName)
                except:
                    badCount += 1
                    badNames.append(tName)
                    continue
                if items[colNums[2]] != 'No hits':
                    dist = float(items[colNums[4]])
                    esd = float(items[colNums[5]])
                else:
                    dist = float(items[colNums[3]])
                    esd = 0.02
                newBond = [[Ids[oInd],Ids[tInd]],['1','1'],dist,esd]
                if newBond not in bondRestData['Bonds']:
                    bondRestData['Bonds'].append(newBond)              
        UpdateBondRestr(bondRestData)
        if badNames:
            msg = f'{badCount} restraints were skipped beccause these atom(s) were not found: {" ".join(set(badNames))}'
            wx.GetApp().Yield()
            G2G.G2MessageBox(G2frame,msg,'Missing atoms')
    def AddAngleRestraint(angleRestData):
        Radii = dict(zip(General['AtomTypes'],zip(General['BondRadii'],General['AngleRadii'])))
        Lists = {'A-atom':[],'B-atom':[],'C-atom':[]}
        for listName in ['A-atom','B-atom']:
            dlg = G2G.G2MultiChoiceDialog(G2frame,'Select '+listName+' for angle A-B-C for '+General['Name']                                                                           ,
                    'Select angle restraint '+listName,Names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                for x in sel:
                    if 'all' in Names[x]:
                        allType = Types[x]
                        for name,Type,coords,Id in zip(Names,Types,Coords,Ids):
                            if Type == allType and 'all' not in name:
                                if 'A' in listName:
                                    Lists[listName].append(Type)
                                else:
                                    Lists[listName].append([Id,Type,coords])
                    else:
                        if 'A' in listName:
                            Lists[listName].append(Types[x])
                        else:
                            Lists[listName].append([Ids[x],Types[x],Coords[x],])
            else:
                return
            targAtoms = [[Ids[x+iBeg],Types[x+iBeg],Coords[x+iBeg]] for x in range(len(Names[iBeg:]))]
        if len(Lists['B-atom']):
            value = 109.54
            dlg = G2G.SingleFloatDialog(G2frame,'Angle','Enter restraint angle ',value,[30.,180.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                value = dlg.GetValue()
            dlg.Destroy()

        Factor = angleRestData['Range']
        indices = (-2,-1,0,1,2)
        Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
        VectA = []
        for Oid,Otype,Ocoord in Lists['B-atom']:
            IndBlist = []
            VectB = []
            for Tid,Ttype,Tcoord in targAtoms:
                result = G2spc.GenAtom(Tcoord,SGData,All=False,Move=False)
                BsumR = (Radii[Otype][0]+Radii[Ttype][0])*Factor
                AsumR = (Radii[Otype][1]+Radii[Ttype][1])*Factor
                for Txyz,Top,Tunit,Spn in result:
                    Dx = (Txyz-Ocoord)+Units
                    dx = np.inner(Amat,Dx)
                    dist = ma.masked_less(np.sqrt(np.sum(dx**2,axis=0)),0.5)
                    IndB = ma.nonzero(ma.masked_greater(dist,BsumR))
                    if np.any(IndB):
                        for indb in IndB:
                            for i in range(len(indb)):
                                if str(dx.T[indb][i]) not in IndBlist:
                                    IndBlist.append(str(dx.T[indb][i]))
                                    unit = Units[indb][i]+Tunit
                                if np.any(unit):
                                    Topstr = '%d+%d,%d,%d'%(Top,unit[0],unit[1],unit[2])
                                else:
                                    Topstr = str(Top)
                                Dist = ma.getdata(dist[indb])[i]
                                if (Dist-AsumR) <= 0.:
                                    VectB.append([Oid,'1',Ocoord,Ttype,Tid,Topstr,Tcoord,Dist])
            VectA.append(VectB)
        for Vects in VectA:
            for i,vecta in enumerate(Vects):                    
                for vectb in Vects[:i]:
                    if vecta[3] in Lists['A-atom']:
                        ids = [vecta[4],vecta[0],vectb[4]]
                        ops = [vecta[5],vecta[1],vectb[5]]
                        angle = [ids,ops,value,1.0]
                        if angle not in angleRestData['Angles']:
                            angleRestData['Angles'].append(angle)
        UpdateAngleRestr(angleRestData)                

    def AddAAAngleRestraint(angleRestData):
        macro = getMacroFile('angle')
        if not macro:
            return
        Chains = makeChains(Names,Ids)            
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Angle']['wtFactor'] = float(items[1])
            elif 'S' in items[0]:
                Res = items[1]
                Value = float(items[2])
                Esd = float(items[3])
                Atms = items[4:7]
                pAtms = ['','','']
                for i,atm in enumerate(Atms):
                    if '+' in atm:
                        pAtms[i] = atm.strip('+')
                ids = [0,0,0]
                chains = list(Chains.keys())
                chains.sort()
                for chain in chains:
                    residues = list(Chains[chain].keys())
                    residues.sort()
                    for residue in residues:
                        for ires in [0,1]:
                            if Res != '*':  #works with disordered res
                                for res,name,Id in Chains[chain][residue][ires]:
                                    if Res == res:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                            else:
                                try:
                                    for res,name,Id in Chains[chain][residue][ires]:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                    for res,name,Id in Chains[chain][residue+1][ires]:
                                        try:
                                            ipos = pAtms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                except KeyError:
                                    continue
                            if all(ids):
                                angle = [list(ids),['1','1','1'],Value,Esd]
                                if angle not in angleRestData['Angles']:
                                    angleRestData['Angles'].append(angle)
                                ids = [0,0,0]
            macStr = macro.readline()
        macro.close()
        print(' Found %d angle restraints'%len(angleRestData['Angles']))
        UpdateAngleRestr(angleRestData)                
        
    def AddMogulAngleRestraint(angleRestData):
        mogul,colNums = getMOGULFile()
        for line in mogul:
            items = line.split(',')
            if 'angle' == items[colNums[0]]:
                aName,bName,cName = items[colNums[1]].split()
                aInd = Names.index(aName)
                bInd = Names.index(bName)
                cInd = Names.index(cName)
                if items[colNums[2]] != 'No hits':
                    angle = float(items[colNums[4]])
                    esd = float(items[colNums[5]])
                else:
                    angle = float(items[colNums[3]])
                    esd = 2.00
                newAngle = [[Ids[aInd],Ids[bInd],Ids[cInd]],['1','1','1'],angle,esd]
                if newAngle not in angleRestData['Angles']:
                    angleRestData['Angles'].append(newAngle)              
        UpdateAngleRestr(angleRestData)                

    def AddPlaneRestraint(restrData):
        ids = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select 4 or more atoms for plane in '+General['Name'],
                'Select 4+ atoms',Names[iBeg:])
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            if len(sel) > 3:
                for x in sel:
                    ids.append(Ids[x+iBeg])
                ops = ['1' for i in range(len(sel))]
                plane = [ids,ops,0.0,0.01]
                if plane not in restrData['Planes']:
                    restrData['Planes'].append(plane)
            else:
                print ('**** ERROR - not enough atoms for a plane restraint - try again ****')
        UpdatePlaneRestr(restrData)                

    def AddAAPlaneRestraint(planeRestData):
        macro = getMacroFile('plane')
        if not macro:
            return
        Chains = makeChains(Names,Ids)            
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Plane']['wtFactor'] = float(items[1])
            elif 'S' in items[0]:
                Res = items[1]
                Esd = float(items[2])
                Atms = items[3:]
                pAtms = ['' for i in Atms]
                for i,atm in enumerate(Atms):
                    if '+' in atm:
                        pAtms[i] = atm.strip('+')
                ids = [0,]*len(Atms)
                ops = ['1' for i in range(len(Atms))]
                chains = list(Chains.keys())
                chains.sort()
                for chain in chains:
                    residues = list(Chains[chain].keys())
                    residues.sort()
                    for residue in residues:
                        for ires in [0,1]:
                            if residue == residues[-1] and Res == '*':
                                continue
                            if Res != '*':  #works with disordered res
                                for res,name,Id in Chains[chain][residue][ires]:
                                    if Res == res:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                            else:
                                try:
                                    for res,name,Id in Chains[chain][residue][ires]:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                    for res,name,Id in Chains[chain][residue+1][ires]:
                                        try:
                                            ipos = pAtms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                except KeyError:
                                    continue
                            if all(ids):
                                plane = [list(ids),ops,0.0,Esd]
                                if plane not in planeRestData['Planes']:
                                    planeRestData['Planes'].append(plane)
                                ids = [0,]*len(Atms)
            macStr = macro.readline()
        macro.close()
        print(' Found %d plane restraints'%len(planeRestData['Planes']))
        UpdatePlaneRestr(planeRestData)                

    def AddAAChiralRestraint(chiralRestData):
        macro = getMacroFile('chiral')
        if not macro:
            return
        Chains = makeChains(Names,Ids)            
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Chiral']['wtFactor'] = float(items[1])
            elif 'S' in items[0]:
                Res = items[1]
                Value = float(items[2])
                Esd = float(items[3])
                Atms = items[4:8]
                ids = [0,0,0,0]
                chains = list(Chains.keys())
                chains.sort()
                for chain in chains:
                    residues = list(Chains[chain].keys())
                    residues.sort()
                    for residue in residues:
                        for ires in [0,1]:
                            if residue == residues[-1] and Res == '*':
                                continue
                            if Res != '*':  #works with disordered res
                                for res,name,Id in Chains[chain][residue][ires]:
                                    if Res == res:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                            else:
                                try:
                                    for res,name,Id in Chains[chain][residue][ires]:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                except KeyError:
                                    continue
                            if all(ids):
                                chiral = [list(ids),['1','1','1','1'],Value,Esd]
                                if chiral not in chiralRestData['Volumes']:
                                    chiralRestData['Volumes'].append(chiral)
                                ids = [0,0,0,0]
            macStr = macro.readline()
        macro.close()
        print(' Found %d chiral volumes restraints'%len(chiralRestData['Volumes']))
        UpdateChiralRestr(chiralRestData)                
        
    def AddAATorsionRestraint(torsionRestData):
        macro = getMacroFile('torsion')
        if not macro:
            return
        Chains = makeChains(Names,Ids)            
        macStr = macro.readline()[:-1]
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Torsion']['wtFactor'] = float(items[1])
            elif 'A' in items[0]:
                name = items[10]
                coeff = np.zeros(9)
                for i,item in enumerate(items[1:10]):
                    coeff[i] = float(item)
                torsionRestData['Coeff'][name] = coeff
            elif 'S' in items[0]:
                Name = items[1]
                Res = items[2]
                Esd = float(items[3])
                Atms = items[4:8]
                pAtms = ['','','','']
                for i,atm in enumerate(Atms):
                    if '+' in atm:
                        pAtms[i] = atm.strip('+')
                ids = [0,0,0,0]
                chains = list(Chains.keys())
                chains.sort()
                for chain in chains:
                    residues = list(Chains[chain].keys())
                    residues.sort()
                    for residue in residues:
                        for ires in [0,1]:
                            if residue == residues[-1] and Res == '*':
                                continue
                            if Res != '*':  #works with disordered res
                                for res,name,Id in Chains[chain][residue][ires]:
                                    if Res == res:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                            else:
                                try:
                                    for res,name,Id in Chains[chain][residue][ires]:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                    for res,name,Id in Chains[chain][residue+1][ires]:
                                        try:
                                            ipos = pAtms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                except KeyError:
                                    continue
                            if all(ids):
                                torsion = [list(ids),['1','1','1','1'],Name,Esd]
                                if torsion not in torsionRestData['Torsions']:
                                    torsionRestData['Torsions'].append(torsion)
                                ids = [0,0,0,0]                            
            macStr = macro.readline()
        macro.close()
        print(' Found %d torsion restraints'%len(torsionRestData['Torsions']))
        UpdateTorsionRestr(torsionRestData)                        
       
    def AddAARamaRestraint(ramaRestData):
        macro = getMacroFile('Ramachandran')
        if not macro:
            return
        Chains = makeChains(Names,Ids)            
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Rama']['wtFactor'] = float(items[1])
            elif 'A' in items[0]:
                nTerms = int(items[1])
                name = items[2]
                coeff = np.zeros((nTerms,6))
                for i in range(nTerms):
                    macStr = macro.readline()
                    items = macStr.split()
                    for j,val in enumerate(items):
                        coeff[i][j] = float(val)
                ramaRestData['Coeff'][name] = coeff
            elif 'S' in items[0]:
                Name = items[1]
                Res = items[2]
                Esd = float(items[3])
                Atms = items[4:9]
                mAtms = ['','','','','']
                pAtms = ['','','','','']
                for i,atm in enumerate(Atms):
                    if '+' in atm:
                        pAtms[i] = atm.strip('+')
                    elif '-' in atm:
                        mAtms[i] = atm.strip('-')
                ids = [0,0,0,0,0]
                chains = list(Chains.keys())
                chains.sort()
                for chain in chains:
                    residues = list(Chains[chain].keys())
                    residues.sort()
                    if not (any(mAtms) or any(pAtms)):
                        for residue in residues:
                            for ires in [0,1]:
                                for res,name,Id in Chains[chain][residue][ires]:
                                    if Res == res:
                                        try:
                                            ipos = Atms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                if all(ids):
                                    rama = [list(ids),['1','1','1','1','1'],Name,Esd]
                                    if rama not in ramaRestData['Ramas']:
                                        ramaRestData['Ramas'].append(rama)
                                    ids = [0,0,0,0,0]
                    else:
                        for ires,residue in enumerate(residues[1:-1]):
                            for jres in [0,1]:
                                try:
                                    for res,name,Id in Chains[chain][residue-1][jres]:
                                        try:
                                            ipos = mAtms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                    for res,name,Id in Chains[chain][residue+1][jres]:
                                        try:
                                            ipos = pAtms.index(name[:3].strip())
                                            ids[ipos] = Id
                                        except ValueError:
                                            continue
                                    for res,name,Id in Chains[chain][residue][jres]:
                                        if Res == res:
                                            try:
                                                ipos = Atms.index(name[:3].strip())
                                                ids[ipos] = Id
                                            except ValueError:
                                                continue
                                    if all(ids):
                                        rama = [list(ids),['1','1','1','1','1'],Name,Esd]
                                        if rama not in ramaRestData['Ramas']:
                                            ramaRestData['Ramas'].append(rama)
                                        ids = [0,0,0,0,0]
                                except KeyError:
                                    continue
            macStr = macro.readline()
        macro.close()
        print(' Found %d Ramachandran restraints'%len(ramaRestData['Ramas']))
        UpdateRamaRestr(ramaRestData)
        
    def AddChemCompRestraint(chemcompRestData):
        ids = []
        factors = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select atoms for chemical restraint in '+General['Name'],
                'Select atoms',Names)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            for x in sel:
                if 'all' in Names[x]:
                    allType = Types[x]
                    for name,Type,Id in zip(Names,Types,Ids):
                        if Type == allType and 'all' not in name:
                            ids.append(Id)
                            factors.append(1.0)
                else:
                    ids.append(Ids[x])
                    factors.append(1.0)
            dlg.Destroy()
            if len(ids) > 0:
                value = 1.0
                dlg = G2G.SingleFloatDialog(G2frame,'Cell sum','Enter unit cell sum ',value,[-1.e6,1.e6],'%.2f')
                if dlg.ShowModal() == wx.ID_OK:
                    value = dlg.GetValue()                
                    comp = [ids,factors,value,0.01]
                    if comp not in chemcompRestData['Sites']:
                        chemcompRestData['Sites'].append(comp)
                UpdateChemcompRestr(chemcompRestData)
            else:
                print ('**** ERROR - not enough atoms for a composition restraint - try again ****')
                
    def AddMomentRestraint(momentRestData):
        ids = []
        dlg = G2G.G2MultiChoiceDialog(G2frame,'Select 2 or more atoms for average moment restraint in '+General['Name'],
                'Select 2+ atoms',Names[iBeg:])
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            if len(sel) > 1:
                for x in sel:
                    ids.append(Ids[x+iBeg])
                moment = [ids,0.0,0.01]
                if moment not in momentRestData['Moments']:
                    momentRestData['Moments'].append(moment)
            else:
                print ('**** ERROR - not enough atoms for a average moment restraint - try again ****')
        UpdateMomentRestr(momentRestData)                

        
    def AddTextureRestraint(textureRestData):
        dlg = wx.TextEntryDialog(G2frame,'Enter h k l for pole figure restraint','Enter HKL','')
        if dlg.ShowModal() == wx.ID_OK:
            text = dlg.GetValue()
            vals = text.split()
            try:
                hkl = [int(vals[i]) for i in range(3)]
                texture = [hkl,24,0.01,False,1.0]
                if texture not in textureRestData['HKLs']:
                        textureRestData['HKLs'].append(texture)
                UpdateTextureRestr(textureRestData)
            except (ValueError,IndexError):
                print ('**** ERROR - bad input of h k l - try again ****')
        dlg.Destroy()
               
    def WtBox(wind,restData):
        if 'Range' not in restData: restData['Range'] = 1.1     #patch
        
        def OnUseData(event):
            Obj = event.GetEventObject()
            restData['Use'] = Obj.GetValue()

        wtBox = wx.BoxSizer(wx.HORIZONTAL)
        wtBox.Add(wx.StaticText(wind,-1,'Phase '+phaseName+' Restraint weight factor: '),0,WACV)
        wtfactor = G2G.ValidatedTxtCtrl(wind,restData,'wtFactor',nDig=(10,2),typeHint=float)
        wtBox.Add(wtfactor,0,WACV)
        useData = wx.CheckBox(wind,-1,label=' Use?')
        useData.Bind(wx.EVT_CHECKBOX, OnUseData)
        useData.SetValue(restData['Use'])        
        wtBox.Add(useData,0,WACV)
        if 'Bonds' in restData or 'Angles' in restData:
            wtBox.Add(wx.StaticText(wind,-1,'  Search range: '),0,WACV)
            sRange = G2G.ValidatedTxtCtrl(wind,restData,'Range',nDig=(10,2),typeHint=float)
            wtBox.Add(sRange,0,WACV)
            wtBox.Add(wx.StaticText(wind,-1,' x sum(atom radii)'),0,WACV)
        return wtBox
        
    def OnRowSelect(event):
        r,c =  event.GetRow(),event.GetCol()
        Obj = event.GetEventObject()
        if r < 0 and c < 0:
            if Obj.IsSelection():
                Obj.ClearSelection()
            else:
                for row in range(Obj.GetNumberRows()):
                    Obj.SelectRow(row,True)
        elif c < 0:                   #only row clicks
            if event.ControlDown():                    
                if r in Obj.GetSelectedRows():
                    Obj.DeselectRow(r)
                else:
                    Obj.SelectRow(r,True)
            elif event.ShiftDown():
                indxs = Obj.GetSelectedRows()
                Obj.ClearSelection()
                ibeg = 0
                if indxs:
                    ibeg = indxs[-1]
                for row in range(ibeg,r+1):
                    Obj.SelectRow(row,True)
            else:
                Obj.ClearSelection()
                Obj.SelectRow(r,True)
                
                    
    def UpdateBondRestr(bondRestData):
        
        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            try:
                new = float(bondTable.GetValue(r,c))
                if new <= 0.:
                    raise ValueError
                bondRestData['Bonds'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdateBondRestr,bondRestData)                

        def OnChangeValue(event):
            rows = GetSelectedRows(Bonds)
            if not rows:
                return
            Bonds.ClearSelection()
            val = bondList[rows[0]][2]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new value for bond',val,[0.,5.],'%.4f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    bondRestData['Bonds'][r][2] = parm
            dlg.Destroy()
            UpdateBondRestr(bondRestData)                

        def OnChangeEsd(event):
            rows = GetSelectedRows(Bonds)
            if not rows:
                return
            Bonds.ClearSelection()
            val = bondList[rows[0]][3]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new esd for bond',val,[0.,1.],'%.4f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    bondRestData['Bonds'][r][3] = parm
            dlg.Destroy()
            UpdateBondRestr(bondRestData)                
                                
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(Bonds,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            Bonds.ClearSelection()
            rows.sort()
            rows.reverse()
            for row in rows:
                bondList.remove(bondList[row])
            UpdateBondRestr(bondRestData)                

        try:
            if BondRestr.GetSizer(): BondRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(BondRestr,bondRestData),0)

        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        bondList = bondRestData['Bonds']
        if len(bondList) and len(bondList[0]) == 6:   #patch
            bondList = bondRestData['Bonds'] = []
        if len(bondList):
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            Types = [wg.GRID_VALUE_STRING,]+4*[wg.GRID_VALUE_FLOAT+':10,3',]
            if 'macro' in General['Type']:
                colLabels = ['(res) A - (res) B','calc','target','esd','delt/sig']
                for i,[indx,ops,obs,esd] in enumerate(bondList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = ''
                        for atom in atoms:
                            name += '('+atom[1]+atom[0].strip()+atom[2]+') '+atom[3]+' - '
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        calc = G2mth.getRestDist(XYZ,Amat)
                        chisq += bondRestData['wtFactor']*((obs-calc)/esd)**2
                        table.append([name[:-3],calc,obs,esd,(obs-calc)/esd])
                        rowLabels.append(str(i))                
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            else:
                colLabels = ['A+SymOp - B+SymOp','calc','target','esd','delt/sig']
                for i,[indx,ops,obs,esd] in enumerate(bondList):
                    try:
                        names = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestDist(XYZ,Amat)
                        chisq += bondRestData['wtFactor']*((obs-calc)/esd)**2
                        table.append([names[0]+'+('+ops[0]+') - '+names[1]+'+('+ops[1]+')',calc,obs,esd,(obs-calc)/esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del bondList[ibad]
            if len(bondList):
                bondTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
                Bonds = G2G.GSGrid(BondRestr)
                Bonds.SetTable(bondTable, True)
                Bonds.AutoSizeColumns(False)
                for c in [0,1,4]: # format read-only rows
                    attr = wx.grid.GridCellAttr()
                    attr.IncRef()
                    attr.SetBackgroundColour(VERY_LIGHT_GREY)
                    attr.SetReadOnly(True)
                    Bonds.SetColAttr(c, attr)
                Bonds.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
                if 'phoenix' in wx.version():
                    Bonds.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
                else:
                    Bonds.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
                for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
                    G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=True)
                G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
                G2frame.Bind(wx.EVT_MENU, OnChangeValue, id=G2G.wxID_RESRCHANGEVAL)
                G2frame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2G.wxID_RESTCHANGEESD)
                frac = '?'
                if 'chisq' in Rvals:
                    frac = f"{100 * chisq / Rvals['chisq']:.1f}"
                mainSizer.Add(wx.StaticText(BondRestr,wx.ID_ANY,
                    f'Bond restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(bondList):.2f}'
                                        ),0)
                Bonds.SetScrollRate(10,10)
                Bonds.SetMinSize((-1,300))
                mainSizer.Add(Bonds,1,wx.EXPAND,1)
            else:
                mainSizer.Add(wx.StaticText(BondRestr,-1,'No bond distance restraints for this phase'),0,)
        else:
            mainSizer.Add(wx.StaticText(BondRestr,-1,'No bond distance restraints for this phase'),0,)

        wx.CallAfter(G2phsGUI.SetPhaseWindow,BondRestr,mainSizer,Scroll=0)
        
    def UpdateAngleRestr(angleRestData):
        
        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            try:
                new = float(angleTable.GetValue(r,c))
                if new <= 0. or new > 180.:
                    raise ValueError
                angleRestData['Angles'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdateAngleRestr,angleRestData)                
            
        def OnChangeValue(event):
            rows = GetSelectedRows(Angles)
            if not rows:
                return
            Angles.ClearSelection()
            val = angleList[rows[0]][2]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new value for angle',val,[0.,360.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    angleRestData['Angles'][r][2] = parm
            dlg.Destroy()
            UpdateAngleRestr(angleRestData)                

        def OnChangeEsd(event):
            rows = GetSelectedRows(Angles)
            if not rows:
                return
            Angles.ClearSelection()
            val = angleList[rows[0]][3]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new esd for angle',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    angleRestData['Angles'][r][3] = parm
            dlg.Destroy()
            UpdateAngleRestr(angleRestData)                
                                            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(Angles,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                angleList.remove(angleList[row])
            UpdateAngleRestr(angleRestData)                
            
        try:
            if AngleRestr.GetSizer():
                AngleRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(AngleRestr,angleRestData),0)

        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        angleList = angleRestData['Angles']
        if len(angleList):
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            Types = [wg.GRID_VALUE_STRING,]+4*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) A - (res) B - (res) C','calc','target','esd','delt/sig']
                for i,[indx,ops,obs,esd] in enumerate(angleList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = ''
                        for atom in atoms:
                            name += '('+atom[1]+atom[0].strip()+atom[2]+') '+atom[3]+' - '
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        calc = G2mth.getRestAngle(XYZ,Amat)
                        chisq += angleRestData['wtFactor']*((obs-calc)/esd)**2
                        table.append([name[:-3],calc,obs,esd,(obs-calc)/esd])
                        rowLabels.append(str(i))                                
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            else:
                colLabels = ['A+SymOp - B+SymOp - C+SymOp','calc','target','esd','delt/sig']
                for i,[indx,ops,obs,esd] in enumerate(angleList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        name = atoms[0]+'+('+ops[0]+') - '+atoms[1]+'+('+ops[1]+') - '+atoms[2]+ \
                        '+('+ops[2]+')'
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestAngle(XYZ,Amat)
                        chisq += angleRestData['wtFactor']*((obs-calc)/esd)**2
                        table.append([name,calc,obs,esd,(obs-calc)/esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del angleList[ibad]
            if len(angleList):
                angleTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
                Angles = G2G.GSGrid(AngleRestr)
                Angles.SetTable(angleTable, True)
                Angles.AutoSizeColumns(False)
                for c in [0,1,4]: # format readonly columns
                    attr = wx.grid.GridCellAttr()
                    attr.IncRef()
                    attr.SetBackgroundColour(VERY_LIGHT_GREY)
                    attr.SetReadOnly(True)
                    Angles.SetColAttr(c, attr)
                Angles.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
                if 'phoenix' in wx.version():
                    Angles.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
                else:
                    Angles.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
                for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
                    G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=True)
                G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
                G2frame.Bind(wx.EVT_MENU, OnChangeValue, id=G2G.wxID_RESRCHANGEVAL)
                G2frame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2G.wxID_RESTCHANGEESD)
                frac = '?'
                if 'chisq' in Rvals:
                    frac = f"{100 * chisq / Rvals['chisq']:.1f}"
                mainSizer.Add(wx.StaticText(AngleRestr,wx.ID_ANY,
                    f'Angle restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(angleList):.2f}'
                                        ),0)
                Angles.SetScrollRate(10,10)
                Angles.SetMinSize((-1,300))
                mainSizer.Add(Angles,1,wx.EXPAND,1)
            else:
                mainSizer.Add(wx.StaticText(AngleRestr,-1,'No bond angle restraints for this phase'),0,)
        else:
            mainSizer.Add(wx.StaticText(AngleRestr,-1,'No bond angle restraints for this phase'),0,)
        wx.CallAfter(G2phsGUI.SetPhaseWindow,AngleRestr,mainSizer,Scroll=0)
    
    def UpdatePlaneRestr(planeRestData):
        
        items = G2frame.dataWindow.RestraintEdit.GetMenuItems()
        for item in items:
            if item.GetItemLabelText() in ['Change value']:
                item.Enable(False)

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            try:
                new = float(planeTable.GetValue(r,c))
                if new <= 0.:
                    raise ValueError
                planeRestData['Planes'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdatePlaneRestr,planeRestData)                
            
        def OnChangeEsd(event):
            rows = GetSelectedRows(Planes)
            if not rows:
                return
            Planes.ClearSelection()
            val = planeList[rows[0]][3]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new esd for plane',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    planeRestData['Planes'][r][3] = parm
            dlg.Destroy()
            UpdatePlaneRestr(planeRestData)                
                                            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(Planes,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                planeList.remove(planeList[row])
            UpdatePlaneRestr(planeRestData)                
            
        try:
            if PlaneRestr.GetSizer():
                PlaneRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(PlaneRestr,planeRestData),0)

        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        planeList = planeRestData['Planes']
        if len(planeList):
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) atom','calc','target','esd']
                for i,[indx,ops,obs,esd] in enumerate(planeList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = ''
                        for a,atom in enumerate(atoms):
                            name += '('+atom[1]+atom[0].strip()+atom[2]+') '+atom[3]+' - '
                            if (a+1)%3 == 0:
                                name += '\n'
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        calc = G2mth.getRestPlane(XYZ,Amat)
                        chisq += planeRestData['wtFactor']*((calc)/esd)**2
                        table.append([name[:-3],calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            else:                                
                colLabels = ['atom+SymOp','calc','target','esd']
                for i,[indx,ops,obs,esd] in enumerate(planeList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestPlane(XYZ,Amat)
                        chisq += planeRestData['wtFactor']*((calc)/esd)**2
                        name = ''
                        for a,atom in enumerate(atoms):
                            name += atom+'+ ('+ops[a]+'),'
                            if (a+1)%3 == 0:
                                name += '\n'
                        table.append([name[:-1],calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del planeList[ibad]
            if len(planeList):
                planeTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
                Planes = G2G.GSGrid(PlaneRestr)
                Planes.SetTable(planeTable, True)
                Planes.AutoSizeColumns(False)
                Planes.AutoSizeRows(False)
                for r in range(len(planeList)):
                    for c in range(3):
                        Planes.SetReadOnly(r,c,True)
                        Planes.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                Planes.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
                if 'phoenix' in wx.version():
                    Planes.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
                else:
                    Planes.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
                for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESTCHANGEESD):
                    G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=True)
                G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
                G2frame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2G.wxID_RESTCHANGEESD)
                frac = '?'
                if 'chisq' in Rvals:
                    frac = f"{100 * chisq / Rvals['chisq']:.1f}"
                mainSizer.Add(wx.StaticText(PlaneRestr,-1,
                    f'Plane restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(planeList):.2f}'
                                        ),0)
                Planes.SetScrollRate(10,10)
                Planes.SetMinSize((-1,300))
                mainSizer.Add(Planes,1,wx.EXPAND,1)
            else:
                mainSizer.Add(wx.StaticText(PlaneRestr,-1,'No plane restraints for this phase'),0,)
        else:
            mainSizer.Add(wx.StaticText(PlaneRestr,-1,'No plane restraints for this phase'),0,)

        G2phsGUI.SetPhaseWindow(PlaneRestr,mainSizer,Scroll=0)
    
    def UpdateChiralRestr(chiralRestData):

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            try:
                new = float(volumeTable.GetValue(r,c))
                if new <= 0.:
                    raise ValueError
                chiralRestData['Volumes'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdateChiralRestr,chiralRestData)                
            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(Volumes,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                volumeList.remove(volumeList[row])
            UpdateChiralRestr(chiralRestData)                
            
        def OnChangeValue(event):
            rows = GetSelectedRows(Volumes)
            if not rows:
                return
            Volumes.ClearSelection()
            val = volumeList[rows[0]][2]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new value for chiral volume',val,[0.,360.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    chiralRestData['Volumes'][r][2] = parm
            dlg.Destroy()
            UpdateChiralRestr(chiralRestData)                

        def OnChangeEsd(event):
            rows = GetSelectedRows(Volumes)
            if not rows:
                return
            Volumes.ClearSelection()
            val = volumeList[rows[0]][3]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new esd for chiral volume',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    chiralRestData['Volumes'][r][3] = parm
            dlg.Destroy()
            UpdateChiralRestr(chiralRestData)                

        try:
            if ChiralRestr.GetSizer(): ChiralRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(ChiralRestr,chiralRestData),0)

        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        volumeList = chiralRestData['Volumes']
        if len(volumeList):
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            Types = [wg.GRID_VALUE_STRING,]+4*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) O (res) A (res) B (res) C','calc','target','esd','delt/sig']
                for i,[indx,ops,obs,esd] in enumerate(volumeList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = ''
                        for atom in atoms:
                            name += '('+atom[1]+atom[0].strip()+atom[2]+') '+atom[3]+' '
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        calc = G2mth.getRestChiral(XYZ,Amat)
                        chisq += chiralRestData['wtFactor']*((obs-calc)/esd)**2
                        table.append([name,calc,obs,esd,(obs-calc)/esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            else:
                colLabels = ['O+SymOp  A+SymOp  B+SymOp  C+SymOp)','calc','target','esd','delt/sig']
                for i,[indx,ops,obs,esd] in enumerate(volumeList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        name = atoms[0]+'+('+ops[0]+') '+atoms[1]+'+('+ops[1]+') '+atoms[2]+ \
                            '+('+ops[2]+') '+atoms[3]+'+('+ops[3]+')'
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestChiral(XYZ,Amat)
                        chisq += chiralRestData['wtFactor']*((obs-calc)/esd)**2
                        table.append([name,calc,obs,esd,(obs-calc)/esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del volumeList[ibad]
            if len(volumeList):
                volumeTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
                Volumes = G2G.GSGrid(ChiralRestr)
                Volumes.SetTable(volumeTable, True)
                Volumes.AutoSizeColumns(False)
                for r in range(len(volumeList)):
                    for c in [0,1,4]:
                        Volumes.SetReadOnly(r,c,True)
                        Volumes.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                Volumes.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
                if 'phoenix' in wx.version():
                    Volumes.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
                else:
                    Volumes.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
                for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
                    G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=True)
                G2frame.Bind(wx.EVT_MENU, OnChangeValue, id=G2G.wxID_RESRCHANGEVAL)
                G2frame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2G.wxID_RESTCHANGEESD)
                G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
                frac = '?'
                if 'chisq' in Rvals:
                    frac = f"{100 * chisq / Rvals['chisq']:.1f}"
                mainSizer.Add(wx.StaticText(ChiralRestr,-1,
                    f'Chiral volume restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(volumeList):.2f}'
                                        ),0)
                Volumes.SetScrollRate(10,10)
                Volumes.SetMinSize((-1,300))
                mainSizer.Add(Volumes,1,wx.EXPAND,1)
            else:
                mainSizer.Add(wx.StaticText(ChiralRestr,-1,'No chiral volume restraints for this phase'),0,)
        else:
            mainSizer.Add(wx.StaticText(ChiralRestr,-1,'No chiral volume restraints for this phase'),0,)

        G2phsGUI.SetPhaseWindow(ChiralRestr,mainSizer,Scroll=0)
    
    def UpdateTorsionRestr(torsionRestData):

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            try:
                new = float(torsionTable.GetValue(r,c))
                if new <= 0. or new > 5.:
                    raise ValueError
                torsionRestData['Torsions'][r][3] = new     #only esd is editable
            except ValueError:
                pass            
            wx.CallAfter(UpdateTorsionRestr,torsionRestData)                
            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(TorsionRestr.Torsions,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                torsionList.remove(torsionList[row])
            wx.CallAfter(UpdateTorsionRestr,torsionRestData)                
            
        def OnChangeEsd(event):
            rows = GetSelectedRows(TorsionRestr.Torsions)
            if not rows:
                return
            TorsionRestr.Torsions.ClearSelection()
            val = torsionList[rows[0]][4]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new esd for torsion restraints',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    torsionRestData['Torsions'][r][4] = parm
            dlg.Destroy()
            wx.CallAfter(UpdateTorsionRestr,torsionRestData) 

        def coeffSizer():               
            table = []
            rowLabels = []
            Types = 9*[wg.GRID_VALUE_FLOAT+':10,4',]
            colLabels = ['Mag A','Pos A','Width A','Mag B','Pos B','Width B','Mag C','Pos C','Width C']
            for item in coeffDict:
                rowLabels.append(item)
                table.append(coeffDict[item])
            coeffTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Coeff = G2G.GSGrid(TorsionRestr)
            Coeff.SetScrollRate(0,10)
            Coeff.SetTable(coeffTable, True)
            Coeff.AutoSizeColumns(False)
            for r in range(len(coeffDict)):
                for c in range(9):
                    Coeff.SetReadOnly(r,c,True)
                    Coeff.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Coeff.SetMaxSize((-1,200))
            return Coeff
                                            
        try:
            if TorsionRestr.GetSizer(): TorsionRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(TorsionRestr,torsionRestData),0)
        
        coeffDict = torsionRestData['Coeff']
        torsionList = torsionRestData['Torsions']
        if len(coeffDict):
            mainSizer.Add((5,5))
            mainSizer.Add(wx.StaticText(TorsionRestr,-1,'Torsion function coefficients:'),0)
            mainSizer.Add(coeffSizer(),1,wx.EXPAND,1)        
        
        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        if len(torsionList):
            mainSizer.Add(wx.StaticText(TorsionRestr,-1,'Torsion restraints:'),0)
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            Types = 2*[wg.GRID_VALUE_STRING,]+4*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) A  B  C  D','coef name','torsion','obs E','restr','esd']
                for i,[indx,ops,cofName,esd] in enumerate(torsionList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = '('+atoms[2][1]+atoms[2][0].strip()+atoms[2][2]+')'
                        for atom in atoms:
                            name += '  '+atom[3]
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        tor = G2mth.getRestTorsion(XYZ,Amat)
                        restr,calc = G2mth.calcTorsionEnergy(tor,coeffDict[cofName])
                        chisq += torsionRestData['wtFactor']*(restr/esd)**2
                        table.append([name,cofName,tor,calc,restr,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del torsionList[ibad]
            torsionTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            TorsionRestr.Torsions = G2G.GSGrid(TorsionRestr)
            TorsionRestr.Torsions.SetTable(torsionTable, True)
            TorsionRestr.Torsions.AutoSizeColumns(False)
            for r in range(len(torsionList)):
                for c in range(5):
                    TorsionRestr.Torsions.SetReadOnly(r,c,True)
                    TorsionRestr.Torsions.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            TorsionRestr.Torsions.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            if 'phoenix' in wx.version():
                TorsionRestr.Torsions.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
            else:
                TorsionRestr.Torsions.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESTCHANGEESD):
                G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=True)
            G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
            G2frame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2G.wxID_RESTCHANGEESD)
            frac = '?'
            if 'chisq' in Rvals:
                frac = f"{100 * chisq / Rvals['chisq']:.1f}"
            mainSizer.Add(wx.StaticText(TorsionRestr,-1,
                f'Torsion restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(torsionList):.2f}'
                                        ),0)
            TorsionRestr.Torsions.SetScrollRate(10,10)
            TorsionRestr.Torsions.SetMinSize((-1,300))
            mainSizer.Add(TorsionRestr.Torsions,1,wx.EXPAND,1)
                
        else:
            mainSizer.Add(wx.StaticText(TorsionRestr,-1,'No torsion restraints for this phase'),0,)

        G2phsGUI.SetPhaseWindow(TorsionRestr,mainSizer,Scroll=0)

    def UpdateRamaRestr(ramaRestData):

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            try:
                new = float(ramaTable.GetValue(r,c))
                if new <= 0. or new > 5.:
                    raise ValueError
                ramaRestData['Ramas'][r][4] = new     #only esd is editable
            except ValueError:
                pass            
            wx.CallAfter(UpdateRamaRestr,ramaRestData)                
            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(RamaRestr.Ramas,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                ramaList.remove(ramaList[row])
            UpdateRamaRestr(ramaRestData)                
            
        def OnChangeEsd(event):
            rows = GetSelectedRows(RamaRestr.Ramas)
            if not rows:
                return
            RamaRestr.Ramas.ClearSelection()
            val = ramaList[rows[0]][4]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new esd for energy',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    ramaRestData['Ramas'][r][4] = parm
            dlg.Destroy()
            UpdateRamaRestr(ramaRestData)

        def coeffSizer():
            table = []
            rowLabels = []
            Types = 6*[wg.GRID_VALUE_FLOAT+':10,4',]
            colLabels = ['Mag','Pos phi','Pos psi','sig(phi)','sig(psi)','sig(cov)']
            for item in coeffDict:
                for i,term in enumerate(coeffDict[item]):
                    rowLabels.append(item+' term:'+str(i))
                    table.append(term)
            coeffTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Coeff = G2G.GSGrid(RamaRestr)
            Coeff.SetScrollRate(0,10)
            Coeff.SetTable(coeffTable, True)
            Coeff.AutoSizeColumns(False)
            for r in range(Coeff.GetNumberRows()):
                for c in range(6):
                    Coeff.SetReadOnly(r,c,True)
                    Coeff.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Coeff.SetMaxSize((-1,200))
            return Coeff
                                                    
        try:
            if RamaRestr.GetSizer(): RamaRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(RamaRestr,ramaRestData),0)
        ramaList = ramaRestData['Ramas']
        coeffDict = ramaRestData['Coeff']
        if len(coeffDict):
            mainSizer.Add(wx.StaticText(RamaRestr,-1,'Ramachandran function coefficients:'),0)
            mainSizer.Add(coeffSizer(),1,wx.EXPAND,1)
            
        if len(ramaList):
            mainSizer.Add(wx.StaticText(RamaRestr,-1,'Ramachandran restraints:'),0)
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            Types = 2*[wg.GRID_VALUE_STRING,]+5*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) A  B  C  D  E','coef name','phi','psi','obs E','restr','esd']
                for i,[indx,ops,cofName,esd] in enumerate(ramaList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = '('+atoms[3][1]+atoms[3][0].strip()+atoms[3][2]+')'
                        for atom in atoms:
                            name += '  '+atom[3]
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        phi,psi = G2mth.getRestRama(XYZ,Amat)
                        restr,calc = G2mth.calcRamaEnergy(phi,psi,coeffDict[cofName])
                        chisq += ramaRestData['wtFactor']*(restr/esd)**2
                        table.append([name,cofName,phi,psi,calc,restr,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print ('**** WARNING - missing atom - restraint deleted ****')
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del ramaList[ibad]
            for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
                G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
            if len(ramaList):
                ramaTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
                RamaRestr.Ramas = G2G.GSGrid(RamaRestr)
                RamaRestr.Ramas.SetTable(ramaTable, True)
                RamaRestr.Ramas.AutoSizeColumns(False)
                for r in range(len(ramaList)):
                    for c in range(6):
                        RamaRestr.Ramas.SetReadOnly(r,c,True)
                        RamaRestr.Ramas.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                RamaRestr.Ramas.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
                if 'phoenix' in wx.version():
                    RamaRestr.Ramas.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
                else:
                    RamaRestr.Ramas.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
                for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESTCHANGEESD):
                    G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=True)
                G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
                G2frame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2G.wxID_RESTCHANGEESD)
                frac = '?'
                if 'chisq' in Rvals:
                    frac = f"{100 * chisq / Rvals['chisq']:.1f}"
                mainSizer.Add(wx.StaticText(RamaRestr,-1,
                    f'Ramachandran restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(ramaList):.2f}'
                                                ),0)
                RamaRestr.Ramas.SetScrollRate(10,10)
                RamaRestr.Ramas.SetMinSize((-1,300))
                mainSizer.Add(RamaRestr.Ramas,1,wx.EXPAND,1)
        else:
            mainSizer.Add(wx.StaticText(RamaRestr,-1,'No Ramachandran restraints for this phase'),0,)

        G2phsGUI.SetPhaseWindow(RamaRestr,mainSizer,Scroll=0)

    def UpdateChemcompRestr(chemcompRestData):
        
        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            rowLabl = ChemComps.GetRowLabelValue(r)
            row = int(rowLabl.split(':')[1])
            if 'Restr' in rowLabl:
                try:
                    new = float(chemcompTable.GetValue(r,c))
                    chemcompRestData['Sites'][row][c-2] = new         #obsd or esd
                except ValueError:
                    pass
            else:
                try:
                    new = float(chemcompTable.GetValue(r,c))
                    id = int(rowLabl.split(':')[2])
                    chemcompRestData['Sites'][row][1][id] = new     #only factor
                except ValueError:
                    pass
            wx.CallAfter(UpdateChemcompRestr,chemcompRestData)                
            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(ChemComps,'delete',G2frame)
            #rows = ChemComps.GetSelectedRows()
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            terms = []
            restrs = []
            for r in sorted(rows,reverse=True):
                rowLabl = ChemComps.GetRowLabelValue(r)                
                if 'Restr' in rowLabl:
                    restrs.append(int(rowLabl.split(':')[1]))
                else:
                    terms.append([int(i) for i in rowLabl.split(':')[1:]])
            # delete terms first in case someone deletes a term in a restraint to be deleted
            for row,term in terms:            
                del chemcompList[row][0][term]
                del chemcompList[row][1][term]
            for row in restrs:
                del chemcompList[row]
            UpdateChemcompRestr(chemcompRestData)                
            
        # def OnChangeValue(event):
        #     rows = GetSelectedRows(ChemComps)
        #     if not rows:
        #         return
        #     ChemComps.ClearSelection()
        #     dlg = G2G.SingleFloatDialog(G2frame,'New value',
        #         'Enter new value for restraint multiplier',1.0,[-1.e6,1.e6],'%.2f')
        #     if dlg.ShowModal() == wx.ID_OK:
        #         parm = dlg.GetValue()
        #         for r in rows:
        #             rowLabl = ChemComps.GetRowLabelValue(r)
        #             if 'term' in rowLabl:
        #                 items = rowLabl.split(':')
        #                 chemcompRestData['Sites'][int(items[1])][1][int(items[2])] = parm
        #     dlg.Destroy()
        #     UpdateChemcompRestr(chemcompRestData)                

        try:
            if ChemCompRestr.GetSizer(): ChemCompRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(ChemCompRestr,chemcompRestData),0)
        mainSizer.Add(wx.StaticText(ChemCompRestr,-1, 
            'NB: The chemical restraint sum is over the unit cell contents'),0)
        mainSizer.Add((5,5),0)

        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        chemcompList = chemcompRestData['Sites']
        if len(chemcompList):
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            Types = [wg.GRID_VALUE_STRING,]+5*[wg.GRID_VALUE_FLOAT+':10,2',]
            colLabels = ['Atoms','mul*frac','factor','calc','target','esd']
            for i,[indx,factors,obs,esd] in enumerate(chemcompList):
                try:
                    atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                    mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs+1))
                    frac = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs-1))
                    mulfrac = mul*frac
                    calcs = mul*frac*factors
                    chisq += chemcompRestData['wtFactor']*((obs-np.sum(calcs))/esd)**2
                    for iatm,[atom,mf,fr,clc] in enumerate(zip(atoms,mulfrac,factors,calcs)):
                        table.append([atom,mf,fr,clc,'',''])
                        rowLabels.append('term:'+str(i)+':'+str(iatm))
                    table.append(['(Sum)','','',np.sum(calcs),obs,esd])
                    rowLabels.append('Restr:'+str(i))
                except KeyError:
                    print ('**** WARNING - missing atom - restraint deleted ****')
                    bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del chemcompList[ibad]
            if len(chemcompList):
                chemcompTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
                ChemComps = G2G.GSGrid(ChemCompRestr)
                ChemComps.SetTable(chemcompTable, True)
                ChemComps.AutoSizeColumns(False)
                for r in range(chemcompTable.GetNumberRows()):
                    for c in range(2):
                        ChemComps.SetReadOnly(r,c,True)
                        ChemComps.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                    if 'Restr' in ChemComps.GetRowLabelValue(r):
                        for c in range(4):
                            ChemComps.SetReadOnly(r,c,True)
                            ChemComps.SetCellStyle(r,c,VERY_YELLOW,True)
                        ChemComps.SetCellTextColour(r,1,VERY_YELLOW) # make spurious #'s disappear
                        ChemComps.SetCellTextColour(r,2,VERY_YELLOW)
                    else:
                        for c in [3,4,5]:
                            ChemComps.SetReadOnly(r,c,True)
                            ChemComps.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                        for c in [4,5]:
                            ChemComps.SetCellTextColour(r,c,VERY_LIGHT_GREY)
                            
                ChemComps.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
                if 'phoenix' in wx.version():
                    ChemComps.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
                else:
                    ChemComps.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
                G2frame.dataWindow.RestraintEdit.Enable(id=G2G.wxID_RESTDELETE,enable=True)
                G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
                #G2frame.Bind(wx.EVT_MENU, OnChangeValue, id=G2G.wxID_RESRCHANGEVAL)
                frac = '?'
                if 'chisq' in Rvals:
                    frac = f"{100 * chisq / Rvals['chisq']:.1f}"
                mainSizer.Add(wx.StaticText(ChemCompRestr,-1,
                    f'Chemical composition restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(chemcompList):.2f}'
                                        ),0)
                mainSizer.Add(ChemComps)
            else:
                mainSizer.Add(wx.StaticText(ChemCompRestr,-1,'No chemical composition restraints for this phase'),0,)
        else:
            mainSizer.Add(wx.StaticText(ChemCompRestr,-1,'No chemical composition restraints for this phase'),0,)

        G2phsGUI.SetPhaseWindow(ChemCompRestr,mainSizer,Scroll=0)
            
    def UpdateMomentRestr(momentRestData):
            
        def OnChangeEsd(event):
            rows = GetSelectedRows(Moments)
            if not rows:
                return
            Moments.ClearSelection()
            val = momentList[rows[0]][-1]
            dlg = G2G.SingleFloatDialog(G2frame,'New value','Enter new esd for moments',val,[0.,5.],'%.3f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    momentRestData['Moments'][r][-1] = parm
            dlg.Destroy()
            UpdateMomentRestr(momentRestData)                
                                            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(Moments,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                momentList.remove(momentList[row])
            UpdateMomentRestr(momentRestData)
                
        try:
            if MomentRestr.GetSizer(): MomentRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(MomentRestr,momentRestData),0)
        mainSizer.Add((5,5),0)

        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
            
        momentList = momentRestData['Moments']
        if len(momentList):
            table = []
            rowLabels = []
            bad = []
            chisq = 0.
            maxno = 0
            for item in momentList:
                maxno = max(maxno,len(item[0]))
            colLabels = maxno*['atom','calc']+['target','esd']
            Types = maxno*[wg.GRID_VALUE_STRING,wg.GRID_VALUE_FLOAT+':10,3']+2*[wg.GRID_VALUE_FLOAT+':10,3',]
            for i,[indx,obs,esd] in enumerate(momentList):
                try:
                    sum = 0.
                    atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                    Mom = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx+4,3))
                    line = []
                    for a,atom in enumerate(atoms):
                        calc = G2mth.GetMag(Mom[a],Cell)
                        sum += calc
                        line += [atom,calc]
                    line += (maxno-len(atoms))*['','']
                    obs = sum/len(atoms)
                    line += [obs,esd]
                    for a,atom in enumerate(atoms):
                        chisq += momentRestData['wtFactor']*((obs-calc)/esd)**2
                    table.append(line)
                    rowLabels.append(str(i))
                except KeyError:
                    print ('**** WARNING - missing atom - restraint deleted ****')
                    bad.append(i)
            momentTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Moments = G2G.GSGrid(MomentRestr)
            Moments.SetTable(momentTable, True)
            Moments.AutoSizeColumns(False)
            Moments.AutoSizeRows(False)
            for r in range(len(momentList)):
                for c in range(maxno*2+1):
                    Moments.SetReadOnly(r,c,True)
                    Moments.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                    if not c%2 and not Moments.GetCellValue(r,c):
                        Moments.SetCellTextColour(r,c+1,VERY_LIGHT_GREY)
            Moments.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESTCHANGEESD):
                G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=True)
            G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
            G2frame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2G.wxID_RESTCHANGEESD)
            frac = '?'
            if 'chisq' in Rvals:
                frac = f"{100 * chisq / Rvals['chisq']:.1f}"
            mainSizer.Add(wx.StaticText(MomentRestr,-1,
                f'Moment restraints: sum(wt*(delt/sig)^2) = {chisq:.2f} ({frac}% of total \u03C7\u00b2), mean(wt*(delt/sig)^2) = {chisq/len(momentList):.2f}'
                                        ),0)
            Moments.SetScrollRate(10,10)
            Moments.SetMinSize((-1,300))
            mainSizer.Add(Moments,1,wx.EXPAND,1)
        else:
            mainSizer.Add(wx.StaticText(MomentRestr,-1,'No magnetic moment restraints for this phase'),0,)
        G2phsGUI.SetPhaseWindow(MomentRestr,mainSizer,Scroll=0)

    def UpdateTextureRestr(textureRestData):
            
        def OnDeleteRestraint(event):
            rows = GetSelectedRows(Textures,'delete',G2frame)
            G2frame.GetStatusBar().SetStatusText('',1)
            if not rows:
                G2frame.GetStatusBar().SetStatusText('First select restraints to be deleted',1)
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                textureList.remove(textureList[row])
            wx.CallAfter(UpdateTextureRestr,textureRestData)                
            
        def OnCellChange(event):
            r,c = event.GetRow(),event.GetCol()
            try:
                if c == 1:  #grid size
                    new = int(textureTable.GetValue(r,c))
                    if new < 6 or new > 24:
                        raise ValueError
                elif c in [2,4]:   #esds
                    new = float(textureTable.GetValue(r,c))
                    if new < -1. or new > 2.:
                        raise ValueError
                else:
                    new = textureTable.GetValue(r,c)
                textureRestData['HKLs'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdateTextureRestr,textureRestData)                

        try:
            if TextureRestr.GetSizer(): TextureRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(TextureRestr,textureRestData),0)
        mainSizer.Add(wx.StaticText(TextureRestr,-1, 
            'NB: The texture restraints suppress negative pole figure values for the selected HKLs\n'
            '    "unit esd" gives a bias toward a flatter polefigure'),0)
        mainSizer.Add((5,5),0)

        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        textureList = textureRestData['HKLs']
        if len(textureList):
            table = []
            rowLabels = []
            Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_LONG,wg.GRID_VALUE_FLOAT+':10,2',
                wg.GRID_VALUE_BOOL,wg.GRID_VALUE_FLOAT+':10,2']
            colLabels = ['HKL','grid','neg esd','use unit?','unit esd']
            for i,[hkl,grid,esd1,ifesd2,esd2] in enumerate(textureList):
                table.append(['%d %d %d'%(hkl[0],hkl[1],hkl[2]),grid,esd1,ifesd2,esd2])
                rowLabels.append(str(i))
            textureTable = G2G.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Textures = G2G.GSGrid(TextureRestr)
            Textures.SetTable(textureTable, True)
            Textures.AutoSizeColumns(False)
            for r in range(len(textureList)):
                Textures.SetReadOnly(r,0,True)
                Textures.SetCellStyle(r,0,VERY_LIGHT_GREY,True)
                if not textureTable.GetValue(r,3):
                    Textures.SetReadOnly(r,4,True)
                    Textures.SetCellStyle(r,4,VERY_LIGHT_GREY,True)
                    Textures.SetCellTextColour(r,4,VERY_LIGHT_GREY)
            Textures.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            if 'phoenix' in wx.version():
                Textures.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
            else:
                Textures.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataWindow.RestraintEdit.Enable(id=G2G.wxID_RESTDELETE,enable=True)
            G2frame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2G.wxID_RESTDELETE)
            mainSizer.Add(Textures,0,)
        else:
            mainSizer.Add(wx.StaticText(TextureRestr,-1,'No texture restraints for this phase'),0,)
        G2phsGUI.SetPhaseWindow(TextureRestr,mainSizer,Scroll=0)
        

    def UpdateGeneralRestr(generalRestData):
        '''Display any generalized restraint expressions'''
        
        def OnEditGenRestraint(event):
            '''Edit a restraint expression'''
            n = event.GetEventObject().index
            parmDict = SetupParmDict(G2frame)
            dlg = G2exG.ExpressionDialog(G2frame,parmDict,
                exprObj=generalRestData['General'][n][0],
                header="Edit a restraint expression",
                fit=False,wildCard=G2frame.testSeqRefineMode())
            restobj = dlg.Show(True)
            if restobj:
                generalRestData['General'][n][0] = restobj
                wx.CallAfter(UpdateGeneralRestr,restrData['General'])
                
        def OnDelGenRestraint(event):              #does this work??
            '''Delete a restraint expression'''
            n = event.GetEventObject().index
            del restrData['General']['General'][n]
            wx.CallAfter(UpdateGeneralRestr,restrData['General'])

        try:
            if GeneralRestr.GetSizer(): GeneralRestr.GetSizer().Clear(True)
        except:
            pass
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        hSizer.Add(wx.StaticText(GeneralRestr,wx.ID_ANY,'Weight factor: '))
        hSizer.Add(G2G.ValidatedTxtCtrl(GeneralRestr,generalRestData,
            'wtFactor',nDig=(10,1),typeHint=float))
        btn = G2G.G2CheckBox(GeneralRestr,'Use?',generalRestData,'Use')
        hSizer.Add(btn)
        hSizer.Add((5,5),0)
        btn = wx.Button(GeneralRestr, wx.ID_ANY,"Add restraint")
        btn.Bind(wx.EVT_BUTTON,OnAddRestraint)
        hSizer.Add(btn,0,wx.EXPAND|wx.ALL)
        mainSizer.Add(hSizer,0)
        mainSizer.Add((5,5),0)
        for i in (G2G.wxID_RESTDELETE,G2G.wxID_RESRCHANGEVAL,G2G.wxID_RESTCHANGEESD):
            G2frame.dataWindow.RestraintEdit.Enable(id=i,enable=False)
        if generalRestData['General']:
            parmDict = SetupParmDict(G2frame)
            GridSiz = wx.FlexGridSizer(0,9,10,2)
            GridSiz.Add((-1,-1))
            GridSiz.Add(
                    wx.StaticText(GeneralRestr,wx.ID_ANY,'Expression'),
                    0,wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL)
            GridSiz.Add((-1,-1))
#            for lbl in ('expression',' ','target\nvalue','current\nvalue','esd'):
            for lbl in ('target\nvalue','current\nvalue','esd'):
                GridSiz.Add(
                    wx.StaticText(GeneralRestr,wx.ID_ANY,lbl,style=wx.CENTER),
                    0,wx.ALIGN_CENTER)
            GridSiz.Add((-1,-1))
            GridSiz.Add((-1,-1))
            GridSiz.Add(
                    wx.StaticText(GeneralRestr,wx.ID_ANY,'Variables',style=wx.CENTER),
                    0,wx.ALIGN_CENTER)
            for i,rest in enumerate(generalRestData['General']):
                eq = rest[0]
                txt = '{}: '.format(i+1)
                GridSiz.Add(wx.StaticText(GeneralRestr,wx.ID_ANY,txt))
                txt = eq.expression
                if len(txt) > 50:
                    txt = txt[:47]+'... '
                txtC = wx.StaticText(GeneralRestr,wx.ID_ANY,txt)
                GridSiz.Add(txtC)
                GridSiz.Add(wx.StaticText(GeneralRestr,wx.ID_ANY,' = '))
                GridSiz.Add(
                    G2G.ValidatedTxtCtrl(GeneralRestr,rest,1,nDig=(10,3,'g'),typeHint=float)
                    )
                # evaluate the expression
                try:
                    calcobj = G2obj.ExpressionCalcObj(rest[0])
                    calcobj.SetupCalc(parmDict)
                    txt = ' {:f} '.format(calcobj.EvalExpression())
                except:
                    txt = ' (error) '
                    txtC.SetForegroundColour("red")
                GridSiz.Add(wx.StaticText(GeneralRestr,wx.ID_ANY,txt))
                GridSiz.Add(
                    G2G.ValidatedTxtCtrl(GeneralRestr,rest,2,nDig=(10,3,'g'),typeHint=float)
                    )
                btn = wx.Button(GeneralRestr, wx.ID_ANY,"Edit",size=(40,-1))
                btn.index = i
                btn.Bind(wx.EVT_BUTTON,OnEditGenRestraint)
                GridSiz.Add(btn)
                btn = wx.Button(GeneralRestr, wx.ID_ANY,"Delete",size=(60,-1))
                btn.index = i
                btn.Bind(wx.EVT_BUTTON,OnDelGenRestraint)
                GridSiz.Add(btn)
                txt = ''
                for i in eq.assgnVars:
                    if txt: txt += '; '
                    txt += str(i) + '=' + str(eq.assgnVars[i])
                if len(txt) > 50:
                    txt = txt[:47]+'...'
                GridSiz.Add(wx.StaticText(GeneralRestr,wx.ID_ANY,txt))
            mainSizer.Add(GridSiz)
        G2phsGUI.SetPhaseWindow(GeneralRestr,mainSizer,Scroll=0)

    
    def OnPageChanged(event):
        page = event.GetSelection()
        #G2frame.restrBook.SetSize(G2frame.dataWindow.GetClientSize())    #TODO -almost right
        text = G2frame.restrBook.GetPageText(page)
        G2frame.dataWindow.RestraintEdit.SetLabel(G2G.wxID_RESRCHANGEVAL,'Change target value')
        G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_USEMOGUL,False)
        if text == 'Bond':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_USEMOGUL,True)
            bondRestData = restrData['Bond']
            UpdateBondRestr(bondRestData)
        elif text == 'Angle':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_USEMOGUL,True)
            angleRestData = restrData['Angle']
            UpdateAngleRestr(angleRestData)
        elif text == 'Plane':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,False)
            planeRestData = restrData['Plane']
            UpdatePlaneRestr(planeRestData)
        elif text == 'Chiral':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,False)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,True)
            chiralRestData = restrData['Chiral']
            UpdateChiralRestr(chiralRestData)
        elif text == 'Torsion':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,False)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,False)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_AARESTRAINTPLOT,True)
            torsionRestData = restrData['Torsion']
            UpdateTorsionRestr(torsionRestData)
        elif text == 'Ramachandran':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,False)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,False)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_AARESTRAINTPLOT,True)
            ramaRestData = restrData['Rama']
            UpdateRamaRestr(ramaRestData)
            wx.CallAfter(G2plt.PlotRama,G2frame,phaseName,rama,ramaName)
        elif text == 'Chem. comp.':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)
            G2frame.dataWindow.RestraintEdit.SetLabel(G2G.wxID_RESRCHANGEVAL,'Change factor')
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTCHANGEESD,False)
            chemcompRestData = restrData['ChemComp']
            UpdateChemcompRestr(chemcompRestData)
        elif text == 'Texture':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,True)
            textureRestData = restrData['Texture']
            UpdateTextureRestr(textureRestData)
        elif text == 'Moments':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,True)
            momentRestData = restrData['Moments']
            UpdateMomentRestr(momentRestData)
        elif text == 'General':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTRAINTADD,True)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESRCHANGEVAL,False)
            G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_RESTCHANGEESD,False)
            UpdateGeneralRestr(restrData['General'])
        event.Skip()

#    def RaisePage(event):
#        'Respond to a "select tab" menu button'
#        # class PseudoEvent(object):
#        #     def __init__(self,page): self.page = page
#        #     def Skip(self): pass
#        #     def GetSelection(self): return self.page
#        try:
#            i = tabIndex.get(event.GetId())
#            G2frame.restrBook.SetSelection(i)
#            #OnPageChanged(PseudoEvent(i))
#        except ValueError:
#            print('Unexpected event in RaisePage')
#
    def OnSelectPage(event):
        'Called when an item is selected from the Select page menu'
        # lookup the menu item that called us and get its text
        tabname = TabSelectionIdDict.get(event.GetId())
        if not tabname:
            print ('Warning: menu item not in dict! id= %d'%event.GetId())
            return
        # find the matching tab
        for PageNum in range(G2frame.restrBook.GetPageCount()):
            if tabname == G2frame.restrBook.GetPageText(PageNum):
                G2frame.restrBook.SetSelection(PageNum)
                return
        else:
            print ("Warning: tab "+tabname+" was not found")

    #### UpdateRestraints execution starts here
    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    lbl= f"Create/edit restraints for {phaseName!r}"[:60]
    topSizer.Add(wx.StaticText(parent,label=lbl),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)
    covdata = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Covariance'))
    Nvars = 0
    Rvals = {}
    if 'Rvals' in covdata:
        Nvars = len(covdata['varyList'])
        Rvals = covdata['Rvals']
        
    try:
        phasedata = G2frame.GetPhaseData()[phaseName]
    except KeyError:        #delete unknown or previously deleted phases from Restraints
        # deleting phase without refresh can cause wx to crash
        #rId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Restraints')
        #pId = G2gd.GetGPXtreeItemId(G2frame,rId,phaseName)
        #G2frame.GPXtree.Delete(pId)
        print('Unknown phase '+phaseName+' is deleted from Restraints')
        mainSizer =  wx.BoxSizer(wx.VERTICAL)
        G2frame.dataWindow.SetSizer(mainSizer)
        mainSizer.Add(
            wx.StaticText(G2frame.dataWindow,-1,' Phase does not exist'),0)
        return
    restrData = data[phaseName]
    if 'Bond' not in restrData:
        restrData['Bond'] = {'wtFactor':1.0,'Range':1.1,'Bonds':[],'Use':True}
    if 'Angle' not in restrData:
        restrData['Angle'] = {'wtFactor':1.0,'Range':0.85,'Angles':[],'Use':True}
    if 'Plane' not in restrData:
        restrData['Plane'] = {'wtFactor':1.0,'Planes':[],'Use':True}
    if 'Chiral' not in restrData:
        restrData['Chiral'] = {'wtFactor':1.0,'Volumes':[],'Use':True}
    if 'Torsion' not in restrData:
        restrData['Torsion'] = {'wtFactor':1.0,'Coeff':{},'Torsions':[],'Use':True}
    if 'Rama' not in restrData:
        restrData['Rama'] = {'wtFactor':1.0,'Coeff':{},'Ramas':[],'Use':True}
    if 'Texture' not in restrData:
        restrData['Texture'] = {'wtFactor':1.0,'HKLs':[],'Use':True}
    if 'ChemComp' not in restrData:
        restrData['ChemComp'] = {'wtFactor':1.0,'Sites':[],'Use':True}
    if 'General' not in restrData:
        restrData['General'] = {'wtFactor':1.0,'General':[], 'Use':True}
    General = phasedata['General']
    if General['Type'] == 'magnetic' and 'Moments' not in restrData:
        restrData['Moments'] = {'wtFactor':1.0,'Moments':[],'Use':True}
    Cell = General['Cell'][1:7]          #skip flag & volume    
    Amat,Bmat = G2lat.cell2AB(Cell)
    SGData = General['SGData']
    cx,ct,cs,cia = General['AtomPtrs']
    Atoms = phasedata['Atoms']
    AtLookUp = G2mth.FillAtomLookUp(Atoms,cia+8)
    if 'macro' in General['Type']:
        Names = [atom[0]+':'+atom[1]+atom[2]+' '+atom[3].ljust(4) for atom in Atoms]
        Ids = []
        Coords = []
        Types = []
        iBeg = 0
    else:    
        Names = ['all '+ name for name in General['AtomTypes']]
        iBeg = len(Names)
        Types = [name for name in General['AtomTypes']]
        Coords = [ [] for type in Types]
        Ids = [ 0 for type in Types]
        Names += [atom[ct-1] for atom in Atoms]
    Types += [atom[ct] for atom in Atoms]
    Coords += [atom[cx:cx+3] for atom in Atoms]
    Ids += [atom[cia+8] for atom in Atoms]
    rama = G2data.ramachandranDist['All']
    ramaName = 'All'
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.RestraintMenu)    
    G2frame.Bind(wx.EVT_MENU, OnAddRestraint, id=G2G.wxID_RESTRAINTADD)
    G2frame.Bind(wx.EVT_MENU, OnUseMogul, id=G2G.wxID_USEMOGUL)
    G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_USEMOGUL,True)
    if 'macro' in phasedata['General']['Type']:
        G2frame.dataWindow.RestraintEdit.Enable(G2G.wxID_AARESTRAINTADD,True)
        G2frame.Bind(wx.EVT_MENU, OnAddAARestraint, id=G2G.wxID_AARESTRAINTADD)
        G2frame.Bind(wx.EVT_MENU, OnPlotAARestraint, id=G2G.wxID_AARESTRAINTPLOT)
    
    # GUI defined here
    G2frame.restrBook = G2G.GSNoteBook(parent=G2frame.dataWindow)
    mainSizer =  wx.BoxSizer(wx.VERTICAL)
    G2frame.dataWindow.SetSizer(mainSizer)
    mainSizer.Add(G2frame.restrBook,1,wx.ALL|wx.EXPAND,1)
    # clear menu and menu pointers
    Pages = []    

    txt = 'Bond'
    BondRestr = wx.ScrolledWindow(G2frame.restrBook)
    G2frame.restrBook.AddPage(BondRestr,txt)
    Pages.append(txt)

    txt = 'Angle'
    AngleRestr = wx.ScrolledWindow(G2frame.restrBook)
    G2frame.restrBook.AddPage(AngleRestr,txt) 
    Pages.append(txt)
   
    txt = 'Plane'
    PlaneRestr = wx.ScrolledWindow(G2frame.restrBook)
    G2frame.restrBook.AddPage(PlaneRestr,txt)
    Pages.append(txt)

    txt = 'Chiral'
    ChiralRestr = wx.ScrolledWindow(G2frame.restrBook)
    G2frame.restrBook.AddPage(ChiralRestr,txt)
    Pages.append(txt)

    if 'macro' in General['Type']:
        txt = 'Torsion'
        TorsionRestr = wx.ScrolledWindow(G2frame.restrBook)
        G2frame.restrBook.AddPage(TorsionRestr,txt)
        Pages.append(txt)

        txt = 'Ramachandran'
        RamaRestr = wx.ScrolledWindow(G2frame.restrBook)
        G2frame.restrBook.AddPage(RamaRestr,txt)
        Pages.append(txt)

    txt = 'Chem. comp.'
    ChemCompRestr = wx.ScrolledWindow(G2frame.restrBook)
    G2frame.restrBook.AddPage(ChemCompRestr,txt)
    Pages.append(txt)
    
    if 'magnetic' in General['Type']:
        txt = 'Moments'
        MomentRestr = wx.ScrolledWindow(G2frame.restrBook)
        G2frame.restrBook.AddPage(MomentRestr,txt)
        Pages.append(txt)
    
    txt = 'General'
    GeneralRestr = wx.ScrolledWindow(G2frame.restrBook)
    G2frame.restrBook.AddPage(GeneralRestr,txt)
    Pages.append(txt)
    
    if General['SH Texture']['Order']:
        txt = 'Texture'
        TextureRestr = wx.ScrolledWindow(G2frame.restrBook)
        G2frame.restrBook.AddPage(TextureRestr,txt)
        Pages.append(txt)

    UpdateBondRestr(restrData['Bond'])

    G2frame.restrBook.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, OnPageChanged)

    # fill page selection menu
    menu = G2frame.dataWindow.RestraintTab
    for page in Pages:
        if menu.FindItem(page) >= 0: continue # is tab already in menu?
        Id = wx.NewId()
        TabSelectionIdDict[Id] = page
        menu.Append(Id,page,'')
        G2frame.Bind(wx.EVT_MENU, OnSelectPage, id=Id)
