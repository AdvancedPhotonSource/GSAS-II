# -*- coding: utf-8 -*-
#GSASIIrestr - restraint GUI routines
########### SVN repository information ###################
# $Date: 2012-12-05 15:38:26 -0600 (Wed, 05 Dec 2012) $
# $Author: vondreele $
# $Revision: 810 $
# $URL: https://subversion.xor.aps.anl.gov/pyGSAS/trunk/GSASIIrestrGUI.py $
# $Id: GSASIIrestrGUI.py 810 2012-12-05 21:38:26Z vondreele $
########### SVN repository information ###################
import wx
import wx.grid as wg
import time
import numpy as np
import numpy.ma as ma
import os.path
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 810 $")
import GSASIImath as G2mth
import GSASIIlattice as G2lat
import GSASIIphsGUI as G2phG
import GSASIIspc as G2spc
import GSASIIgrid as G2gd
import GSASIIplot as G2plt
import GSASIIdata as G2data

VERY_LIGHT_GREY = wx.Colour(235,235,235)

################################################################################
#####  Restraints
################################################################################           
       
def UpdateRestraints(G2frame,data,Phases,phaseName):
    if not len(Phases):
        print 'There are no phases to form restraints'
        return
    phasedata = Phases[phaseName]
    if phaseName not in data:
        data[phaseName] = {}
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
    General = phasedata['General']
    Cell = General['Cell'][1:7]          #skip flag & volume    
    Amat,Bmat = G2lat.cell2AB(Cell)
    SGData = General['SGData']
    cx,ct,cs = General['AtomPtrs'][:3]
    Atoms = phasedata['Atoms']
    AtLookUp = G2mth.FillAtomLookUp(Atoms)
    if 'macro' in General['Type']:
        Names = [atom[0]+':'+atom[1]+atom[2]+' '+atom[3] for atom in Atoms]
        Ids = []
        Coords = []
        Types = []
    else:    
        Names = ['all '+ name for name in General['AtomTypes']]
        iBeg = len(Names)
        Types = [name for name in General['AtomTypes']]
        Coords = [ [] for type in Types]
        Ids = [ 0 for type in Types]
        Names += [atom[ct-1] for atom in Atoms]
    Types += [atom[ct] for atom in Atoms]
    Coords += [atom[cx:cx+3] for atom in Atoms]
    Ids += [atom[-1] for atom in Atoms]
    rama = G2data.ramachandranDist['All']
    ramaName = 'All'
    
    def OnSelectPhase(event):
        dlg = wx.SingleChoiceDialog(G2frame,'Select','Phase',Phases.keys())
        try:
            if dlg.ShowModal() == wx.ID_OK:
                phaseName = Phases.keys()[dlg.GetSelection()]
                UpdateRestraints(G2frame,data,Phases,phaseName)
        finally:
            dlg.Destroy()
    
    def getMacroFile(macName):
        defDir = os.path.join(os.path.split(__file__)[0],'GSASIImacros')
        dlg = wx.FileDialog(G2frame,message='Choose '+macName+' restraint macro file',
            defaultDir=defDir,defaultFile="",wildcard="GSAS-II macro file (*.mac)|*.mac",
            style=wx.OPEN | wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                macfile = dlg.GetPath()
                macro = open(macfile,'Ur')
                head = macro.readline()
                if macName not in head:
                    print head
                    print '**** ERROR - wrong restraint macro file selected, try again ****'
                    macro = []
            else: # cancel was pressed
                macxro = []
        finally:
            dlg.Destroy()
        return macro        #advanced past 1st line
        
    def OnPlotAARestraint(event):
        page = G2frame.dataDisplay.GetSelection()
        if 'Torsion' in G2frame.dataDisplay.GetPageText(page):
            torNames = []
            torNames += restrData['Torsion']['Coeff'].keys()
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
            
        elif 'Rama' in G2frame.dataDisplay.GetPageText(page):
            ramaNames = ['All',]
            ramaNames += restrData['Rama']['Coeff'].keys()
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
            
    def OnAddRestraint(event):
        page = G2frame.dataDisplay.GetSelection()
        if 'Bond' in G2frame.dataDisplay.GetPageText(page):
            AddBondRestraint(restrData['Bond'])
        elif 'Angle' in G2frame.dataDisplay.GetPageText(page):
            AddAngleRestraint(restrData['Angle'])
        elif 'Plane' in G2frame.dataDisplay.GetPageText(page):
            AddPlaneRestraint(restrData['Plane'])
        elif 'Chem' in G2frame.dataDisplay.GetPageText(page):
            AddChemCompRestraint(restrData['ChemComp'])
        elif 'Texture' in G2frame.dataDisplay.GetPageText(page):
            AddTextureRestraint(restrData['Texture'])
            
    def OnAddAARestraint(event):
        page = G2frame.dataDisplay.GetSelection()
        if 'Bond' in G2frame.dataDisplay.GetPageText(page):
            AddAABondRestraint(restrData['Bond'])
        elif 'Angle' in G2frame.dataDisplay.GetPageText(page):
            AddAAAngleRestraint(restrData['Angle'])
        elif 'Plane' in G2frame.dataDisplay.GetPageText(page):
            AddAAPlaneRestraint(restrData['Plane'])
        elif 'Chiral' in G2frame.dataDisplay.GetPageText(page):
            AddAAChiralRestraint(restrData['Chiral'])
        elif 'Torsion' in G2frame.dataDisplay.GetPageText(page):
            AddAATorsionRestraint(restrData['Torsion'])
        elif 'Rama' in G2frame.dataDisplay.GetPageText(page):
            AddAARamaRestraint(restrData['Rama'])
            
    def AddBondRestraint(bondRestData):
        Lists = {'origin':[],'target':[]}
        for listName in ['origin','target']:
            dlg = wx.MultiChoiceDialog(G2frame,'Bond restraint '+listName+' for '+General['Name'],
                    'Select bond restraint '+listName+' atoms',Names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                for x in sel:
                    if 'all' in Names[x]:
                        allType = Types[x]
                        for name,Type,coords,id in zip(Names,Types,Coords,Ids):
                            if Type == allType and 'all' not in name:
                                Lists[listName].append([id,Type,coords])
                    else:
                        Lists[listName].append([Ids[x],Types[x],Coords[x],])
            else:
                break
        if len(Lists['origin']) and len(Lists['target']):
            bond = 1.54
            dlg = G2phG.SingleFloatDialog(G2frame,'Distance','Enter restraint distance for bond',bond,[0.01,4.],'%.4f')
            if dlg.ShowModal() == wx.ID_OK:
                bond = dlg.GetValue()
            dlg.Destroy()
        Factor = bondRestData['Range']
        indices = (-1,0,1)
        Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
        origAtoms = Lists['origin']
        targAtoms = Lists['target']
        dlg = wx.ProgressDialog("Generating bond restraints","Processed origin atoms",len(origAtoms), 
            style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME)
        try:
            Norig = 0
            for Oid,Otype,Ocoord in origAtoms:
                Norig += 1
                dlg.Update(Norig)
                for Tid,Ttype,Tcoord in targAtoms:
                    if 'macro' in General['Type']:
                        result = [[Tcoord,1,[0,0,0]],]
                    else:
                        result = G2spc.GenAtom(Tcoord,SGData,False,Move=False)
                    for Txyz,Top,Tunit in result:
                        Dx = (Txyz-np.array(Ocoord))+Units
                        dx = np.inner(Amat,Dx)
                        dist = ma.masked_less(np.sqrt(np.sum(dx**2,axis=0)),bond/Factor)
                        IndB = ma.nonzero(ma.masked_greater(dist,bond*Factor))
                        if np.any(IndB):
                            for indb in IndB:
                                for i in range(len(indb)):
                                    unit = Units[indb][i]+Tunit
                                    if np.any(unit):
                                        Topstr = '%d+%d,%d,%d'%(Top,unit[0],unit[1],unit[2])
                                    else:
                                        Topstr = str(Top)
                                    newBond = [[Oid,Tid],['1',Topstr],bond,0.01]
                                    if newBond not in bondRestData['Bonds']:
                                        bondRestData['Bonds'].append(newBond)
        finally:
            dlg.Destroy()
        UpdateBondRestr(bondRestData)                

    def AddAABondRestraint(bondRestData):
        Radii = dict(zip(General['AtomTypes'],General['BondRadii']))
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
                tIds = []
                tCoords = []
                res = items[1]
                dist = float(items[2])
                esd = float(items[3])
                oAtm,tAtm = items[4:6]
                for Name,coords,Id in atoms:
                    names = Name.split()
                    if res == '*' or res in names[0]:
                        if oAtm == names[2]:
                            oIds.append(Id)
                            oCoords.append(np.array(coords))
                        if tAtm == names[2]:
                            tIds.append(Id)
                            tCoords.append(np.array(coords))
                for i,[oId,oCoord] in enumerate(zip(oIds,oCoords)):
                    for tId,tCoord in zip(tIds,tCoords)[i:]:
                        obsd = np.sqrt(np.sum(np.inner(Amat,tCoord-oCoord)**2))
                        if dist/Factor < obsd < dist*Factor:
                            newBond = [[oId,tId],['1','1'],dist,esd]
                            if newBond not in bondRestData['Bonds']:
                                bondRestData['Bonds'].append(newBond)                          
            macStr = macro.readline()
        macro.close()
        UpdateBondRestr(bondRestData)                
            
    def AddAngleRestraint(angleRestData):
        Radii = dict(zip(General['AtomTypes'],zip(General['BondRadii'],General['AngleRadii'])))
        Lists = {'A-atom':[],'B-atom':[],'C-atom':[]}
        for listName in ['A-atom','B-atom']:
            dlg = wx.MultiChoiceDialog(G2frame,'Select '+listName+' for angle A-B-C for '+General['Name']                                                                           ,
                    'Select angle restraint '+listName,Names)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelections()
                for x in sel:
                    if 'all' in Names[x]:
                        allType = Types[x]
                        for name,Type,coords,id in zip(Names,Types,Coords,Ids):
                            if Type == allType and 'all' not in name:
                                if 'A' in listName:
                                    Lists[listName].append(Type)
                                else:
                                    Lists[listName].append([id,Type,coords])
                    else:
                        if 'A' in listName:
                            Lists[listName].append(Types[x])
                        else:
                            Lists[listName].append([Ids[x],Types[x],Coords[x],])
            else:
                break
            targAtoms = [[Ids[x+iBeg],Types[x+iBeg],Coords[x+iBeg]] for x in range(len(Names[iBeg:]))]
        if len(Lists['B-atom']):
            value = 109.54
            dlg = G2phG.SingleFloatDialog(G2frame,'Angle','Enter restraint angle ',value,[30.,180.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                value = dlg.GetValue()
            dlg.Destroy()

        Factor = angleRestData['Range']
        indices = (-1,0,1)
        Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
        VectA = []
        for Oid,Otype,Ocoord in Lists['B-atom']:
            IndBlist = []
            VectB = []
            for Tid,Ttype,Tcoord in targAtoms:
                result = G2spc.GenAtom(Tcoord,SGData,False,Move=False)
                BsumR = (Radii[Otype][0]+Radii[Ttype][0])*Factor
                AsumR = (Radii[Otype][1]+Radii[Ttype][1])*Factor
                for Txyz,Top,Tunit in result:
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
                        XYZ = np.array([vecta[6],vecta[2],vectb[6]])
                        angle = [ids,ops,value,1.0]
                        if angle not in angleRestData['Angles']:
                            angleRestData['Angles'].append(angle)
        UpdateAngleRestr(angleRestData)                

    def AddAAAngleRestraint(angleRestData):
        macro = getMacroFile('angle')
        if not macro:
            return
        atoms = zip(Names,Ids)
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Angle']['wtFactor'] = float(items[1])
            elif 'S' in items[0]:
                res = items[1]
                value = float(items[2])
                esd = float(items[3])
                Atms = items[4:7]
                pAtms = ['','','']
                for i,atm in enumerate(Atms):
                    if '+' in atm:
                        pAtms[i] = atm.strip('+')
                ids = np.array([0,0,0])
                rNum = -1
                for name,id in atoms:
                    names = name.split()
                    tNum = int(names[0].split(':')[0])
                    if res in names[0]:
                        try:
                            ipos = Atms.index(names[2])
                            ids[ipos] = id
                        except ValueError:
                            continue
                    elif res == '*':
                        try:
                            ipos = Atms.index(names[2])
                            if not np.all(ids):
                                rNum = int(names[0].split(':')[0])
                            ids[ipos] = id
                        except ValueError:
                            try:
                                if tNum == rNum+1:
                                    ipos = pAtms.index(names[2])
                                    ids[ipos] = id
                            except ValueError:
                                continue
                    if np.all(ids):
                        angle = [list(ids),['1','1','1'],value,esd]
                        if angle not in angleRestData['Angles']:
                            angleRestData['Angles'].append(angle)
                        ids = np.array([0,0,0])
            macStr = macro.readline()
        macro.close()
        UpdateAngleRestr(angleRestData)                
        
    def AddPlaneRestraint(restrData):
        ids = []
        dlg = wx.MultiChoiceDialog(G2frame,'Select 4 or more atoms for plane in '+General['Name'],
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
                print '**** ERROR - not enough atoms for a plane restraint - try again ****'
        UpdatePlaneRestr(restrData)                

    def AddAAPlaneRestraint(planeRestData):
        macro = getMacroFile('plane')
        if not macro:
            return
        atoms = zip(Names,Ids)
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Plane']['wtFactor'] = float(items[1])
            elif 'S' in items[0]:
                res = items[1]
                esd = float(items[2])
                Atms = items[3:]
                pAtms = ['' for i in Atms]
                for i,atm in enumerate(Atms):
                    if '+' in atm:
                        pAtms[i] = atm.strip('+')
                rNum = -1
                ids = np.zeros(len(Atms))
                ops = ['1' for i in range(len(Atms))]
                for name,id in atoms:
                    names = name.split()
                    tNum = int(names[0].split(':')[0])
                    if res in names[0]:
                        try:
                            ipos = Atms.index(names[2])
                            ids[ipos] = id
                        except ValueError:
                            continue
                    elif res == '*':
                        try:
                            ipos = Atms.index(names[2])
                            if not np.all(ids):
                                rNum = int(names[0].split(':')[0])
                            ids[ipos] = id
                        except ValueError:
                            try:
                                if tNum == rNum+1:
                                    ipos = pAtms.index(names[2])
                                    ids[ipos] = id
                            except ValueError:
                                continue
                    if np.all(ids):
                        plane = [list(ids),ops,0.0,esd]
                        if plane not in planeRestData['Planes']:
                            planeRestData['Planes'].append(plane)
                        ids = np.zeros(len(Atms))
            macStr = macro.readline()
        macro.close()
        UpdatePlaneRestr(planeRestData)                

    def AddAAChiralRestraint(chiralRestData):
        macro = getMacroFile('chiral')
        if not macro:
            return
        atoms = zip(Names,Ids)
        macStr = macro.readline()
        while macStr:
            items = macStr.split()
            if 'F' in items[0]:
                restrData['Chiral']['wtFactor'] = float(items[1])
            elif 'S' in items[0]:
                res = items[1]
                value = float(items[2])
                esd = float(items[3])
                Atms = items[4:8]
                ids = np.array([0,0,0,0])
                for name,id in atoms:
                    names = name.split()
                    if res in names[0]:
                        try:
                            ipos = Atms.index(names[2])
                            ids[ipos] = id
                        except ValueError:
                            pass
                        if np.all(ids):
                            chiral = [list(ids),['1','1','1','1'],value,esd]
                            if chiral not in chiralRestData['Volumes']:
                                chiralRestData['Volumes'].append(chiral)
                            ids = np.array([0,0,0,0])
            macStr = macro.readline()
        macro.close()
        UpdateChiralRestr(chiralRestData)                
        
    def makeChains(Names,Ids):
        Chains = {}
        chain = ''
        atoms = zip(Names,Ids)
        for name,id in atoms:
            items = name.split()
            rnum,res = items[0].split(':')
            if items[1] not in Chains:
                Residues = {}
                Chains[items[1]] = Residues
            if int(rnum) not in Residues:
                Residues[int(rnum)] = []
            Residues[int(rnum)].append([res,items[2],id])
        return Chains
        
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
                ids = np.array([0,0,0,0])
                chains = Chains.keys()
                chains.sort()
                for chain in chains:
                    residues = Chains[chain].keys()
                    residues.sort()
                    for residue in residues:
                        if residue == residues[-1] and Res == '*':
                            continue
                        if Res != '*':
                            for res,name,id in Chains[chain][residue]:
                                if Res == res:
                                    try:
                                        ipos = Atms.index(name)
                                        ids[ipos] = id
                                    except ValueError:
                                        continue
                        else:
                            for res,name,id in Chains[chain][residue]:
                                try:
                                    ipos = Atms.index(name)
                                    ids[ipos] = id
                                except ValueError:
                                    continue
                            for res,name,id in Chains[chain][residue+1]:
                                try:
                                    ipos = pAtms.index(name)
                                    ids[ipos] = id
                                except ValueError:
                                    continue
                        if np.all(ids):
                            torsion = [list(ids),['1','1','1','1'],Name,Esd]
                            if torsion not in torsionRestData['Torsions']:
                                torsionRestData['Torsions'].append(torsion)
                            ids = np.array([0,0,0,0])                            
            macStr = macro.readline()
        macro.close()
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
                ids = np.array([0,0,0,0,0])
                chains = Chains.keys()
                chains.sort()
                for chain in chains:
                    residues = Chains[chain].keys()
                    residues.sort()
                    if not (any(mAtms) or any(pAtms)):
                        for residue in residues:
                            for res,name,id in Chains[chain][residue]:
                                if Res == res:
                                    try:
                                        ipos = Atms.index(name)
                                        ids[ipos] = id
                                    except ValueError:
                                        continue
                            if np.all(ids):
                                rama = [list(ids),['1','1','1','1','1'],Name,Esd]
                                if rama not in ramaRestData['Ramas']:
                                    ramaRestData['Ramas'].append(rama)
                                ids = np.array([0,0,0,0,0])
                    else:
                        for residue in residues[1:-1]:
                            for res,name,id in Chains[chain][residue-1]:
                                try:
                                    ipos = mAtms.index(name)
                                    ids[ipos] = id
                                except ValueError:
                                    continue
                            for res,name,id in Chains[chain][residue+1]:
                                try:
                                    ipos = pAtms.index(name)
                                    ids[ipos] = id
                                except ValueError:
                                    continue
                            for res,name,id in Chains[chain][residue]:
                                if Res == res:
                                    try:
                                        ipos = Atms.index(name)
                                        ids[ipos] = id
                                    except ValueError:
                                        continue
                            if np.all(ids):
                                rama = [list(ids),['1','1','1','1','1'],Name,Esd]
                                if rama not in ramaRestData['Ramas']:
                                    ramaRestData['Ramas'].append(rama)
                                ids = np.array([0,0,0,0,0])
            macStr = macro.readline()
        macro.close()
        UpdateRamaRestr(ramaRestData)
        
    def AddChemCompRestraint(chemcompRestData):
        ids = []
        factors = []
        dlg = wx.MultiChoiceDialog(G2frame,'Select atoms for chemical restraint in '+General['Name'],
                'Select atoms',Names)
        if dlg.ShowModal() == wx.ID_OK:
            sel = dlg.GetSelections()
            for x in sel:
                if 'all' in Names[x]:
                    allType = Types[x]
                    for name,Type,id in zip(Names,Types,Ids):
                        if Type == allType and 'all' not in name:
                            ids.append(id)
                            factors.append(1.0)
                else:
                    ids.append(Ids[x])
                    factors.append(1.0)
            dlg.Destroy()
            if len(ids) > 2:
                value = 1.0
                dlg = G2phG.SingleFloatDialog(G2frame,'Angle','Enter unit cell sum ',value,[-1.e6,1.e6],'%.2f')
                if dlg.ShowModal() == wx.ID_OK:
                    value = dlg.GetValue()                
                    comp = [ids,factors,value,0.01]
                    if comp not in chemcompRestData['Sites']:
                        chemcompRestData['Sites'].append(comp)
                UpdateChemcompRestr(chemcompRestData)
            else:
                print '**** ERROR - not enough atoms for a composition restraint - try again ****'
        
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
                print '**** ERROR - bad input of h k l - try again ****'
        dlg.Destroy()
               
    def WtBox(wind,restData):
        if 'Range' not in restData: restData['Range'] = 1.1     #patch
        
        def OnWtFactor(event):
            try:
                value = float(wtfactor.GetValue())
            except ValueError:
                value = 1.0
            restData['wtFactor'] = value
            wtfactor.SetValue('%.2f'%(value))
            
        def OnRange(event):
            try:
                value = float(sRange.GetValue())
            except ValueError:
                value = 1.0
            restData['Range'] = value
            sRange.SetValue('%.2f'%(value))
            
        def OnUseData(event):
            Obj = event.GetEventObject()
            restData['Use'] = Obj.GetValue()

        wtBox = wx.BoxSizer(wx.HORIZONTAL)
        wtBox.Add(wx.StaticText(wind,-1,'Restraint weight factor:'),0,wx.ALIGN_CENTER_VERTICAL)
        wtfactor = wx.TextCtrl(wind,-1,value='%.2f'%(restData['wtFactor']),style=wx.TE_PROCESS_ENTER)
        wtfactor.Bind(wx.EVT_TEXT_ENTER,OnWtFactor)
        wtfactor.Bind(wx.EVT_KILL_FOCUS,OnWtFactor)
        wtBox.Add(wtfactor,0,wx.ALIGN_CENTER_VERTICAL)
        useData = wx.CheckBox(wind,-1,label=' Use?')
        useData.Bind(wx.EVT_CHECKBOX, OnUseData)
        useData.SetValue(restData['Use'])        
        wtBox.Add(useData,0,wx.ALIGN_CENTER_VERTICAL)
        if 'Bonds' in restData or 'Angles' in restData:
            wtBox.Add(wx.StaticText(wind,-1,'Search range:'),0,wx.ALIGN_CENTER_VERTICAL)
            sRange = wx.TextCtrl(wind,-1,value='%.2f'%(restData['Range']),style=wx.TE_PROCESS_ENTER)
            sRange.Bind(wx.EVT_TEXT_ENTER,OnRange)
            sRange.Bind(wx.EVT_KILL_FOCUS,OnRange)
            wtBox.Add(sRange,0,wx.ALIGN_CENTER_VERTICAL)
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
            val = bondRestData['Bonds'][r][c]
            try:
                new = float(bondTable.GetValue(r,c))
                if new <= 0.:
                    raise ValueError
                bondRestData['Bonds'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdateBondRestr,bondRestData)                
            
        def OnChangeValue(event):
            rows = Bonds.GetSelectedRows()
            if not rows:
                return
            Bonds.ClearSelection()
            val = bondList[rows[0]][2]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new value for bond',val,[0.,5.],'%.4f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    bondRestData['Bonds'][r][2] = parm
            dlg.Destroy()
            UpdateBondRestr(bondRestData)                

        def OnChangeEsd(event):
            rows = Bonds.GetSelectedRows()
            if not rows:
                return
            Bonds.ClearSelection()
            val = bondList[rows[0]][3]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for bond',val,[0.,1.],'%.4f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    bondRestData['Bonds'][r][3] = parm
            dlg.Destroy()
            UpdateBondRestr(bondRestData)                
                                
        def OnDeleteRestraint(event):
            rows = Bonds.GetSelectedRows()
            if not rows:
                return
            Bonds.ClearSelection()
            rows.sort()
            rows.reverse()
            for row in rows:
                bondList.remove(bondList[row])
            UpdateBondRestr(bondRestData)                
            
        BondRestr.DestroyChildren()
        dataDisplay = wx.Panel(BondRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(BondRestr,bondRestData),0,wx.ALIGN_CENTER_VERTICAL)

        bondList = bondRestData['Bonds']
        if len(bondList) and len(bondList[0]) == 6:   #patch
            bondList = bondRestData['Bonds'] = []
        if len(bondList):
            table = []
            rowLabels = []
            bad = []
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,3',]
            if 'macro' in General['Type']:
                colLabels = ['(res) A - (res) B','calc','obs','esd']
                for i,[indx,ops,obs,esd] in enumerate(bondList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = ''
                        for atom in atoms:
                            name += '('+atom[1]+atom[0].strip()+atom[2]+') '+atom[3]+' - '
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        calc = G2mth.getRestDist(XYZ,Amat)
                        table.append([name[:-3],calc,obs,esd])
                        rowLabels.append(str(i))                
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            else:
                colLabels = ['A+SymOp - B+SymOp','calc','obs','esd']
                for i,[indx,ops,obs,esd] in enumerate(bondList):
                    try:
                        names = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestDist(XYZ,Amat)
                        table.append([names[0]+'+('+ops[0]+') - '+names[1]+'+('+ops[1]+')',calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del bondList[ibad]
            bondTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Bonds = G2gd.GSGrid(BondRestr)
            Bonds.SetTable(bondTable, True)
            Bonds.AutoSizeColumns(False)
            for r in range(len(bondList)):
                for c in range(2):
                    Bonds.SetReadOnly(r,c,True)
                    Bonds.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Bonds.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            Bonds.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeValue, id=G2gd.wxID_RESRCHANGEVAL)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2gd.wxID_RESTCHANGEESD)
            mainSizer.Add(Bonds,0,)
        else:
            mainSizer.Add(wx.StaticText(BondRestr,-1,'No bond distance restraints for this phase'),0,)

        BondRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50       #make room for tab
        BondRestr.SetSize(Size)
        BondRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)
        
    def UpdateAngleRestr(angleRestData):
        
        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            val = angleRestData['Angles'][r][c]
            try:
                new = float(angleTable.GetValue(r,c))
                if new <= 0. or new > 180.:
                    raise ValueError
                angleRestData['Angles'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdateAngleRestr,angleRestData)                
            
        def OnChangeValue(event):
            rows = Angles.GetSelectedRows()
            if not rows:
                return
            Angles.ClearSelection()
            val = angleList[rows[0]][2]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new value for angle',val,[0.,360.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    angleRestData['Angles'][r][2] = parm
            dlg.Destroy()
            UpdateAngleRestr(angleRestData)                

        def OnChangeEsd(event):
            rows = Angles.GetSelectedRows()
            if not rows:
                return
            Angles.ClearSelection()
            val = angleList[rows[0]][3]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for angle',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    angleRestData['Angles'][r][3] = parm
            dlg.Destroy()
            UpdateAngleRestr(angleRestData)                
                                            
        def OnDeleteRestraint(event):
            rows = Angles.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                angleList.remove(angleList[row])
            UpdateAngleRestr(angleRestData)                
            
        AngleRestr.DestroyChildren()
        dataDisplay = wx.Panel(AngleRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(AngleRestr,angleRestData),0,wx.ALIGN_CENTER_VERTICAL)

        angleList = angleRestData['Angles']
        if len(angleList):
            table = []
            rowLabels = []
            bad = []
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) A - (res) B - (res) C','calc','obs','esd']
                for i,[indx,ops,obs,esd] in enumerate(angleList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = ''
                        for atom in atoms:
                            name += '('+atom[1]+atom[0].strip()+atom[2]+') '+atom[3]+' - '
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        calc = G2mth.getRestAngle(XYZ,Amat)
                        table.append([name[:-3],calc,obs,esd])
                        rowLabels.append(str(i))                                
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            else:
                colLabels = ['A+SymOp - B+SymOp - C+SymOp','calc','obs','esd']
                for i,[indx,ops,obs,esd] in enumerate(angleList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        name = atoms[0]+'+('+ops[0]+') - '+atoms[1]+'+('+ops[1]+') - '+atoms[2]+ \
                        '+('+ops[2]+')'
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestAngle(XYZ,Amat)
                        table.append([name,calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del angleList[ibad]
            angleTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Angles = G2gd.GSGrid(AngleRestr)
            Angles.SetTable(angleTable, True)
            Angles.AutoSizeColumns(False)
            for r in range(len(angleList)):
                for c in range(2):
                    Angles.SetReadOnly(r,c,True)
                    Angles.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Angles.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            Angles.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeValue, id=G2gd.wxID_RESRCHANGEVAL)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2gd.wxID_RESTCHANGEESD)
            mainSizer.Add(Angles,0,)
        else:
            mainSizer.Add(wx.StaticText(AngleRestr,-1,'No bond angle restraints for this phase'),0,)

        AngleRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50      #make room for tab
        AngleRestr.SetSize(Size)
        AngleRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)
    
    def UpdatePlaneRestr(planeRestData):
        
        items = G2frame.dataFrame.RestraintEdit.GetMenuItems()
        for item in items:
            if item.GetLabel() in ['Change value']:
                item.Enable(False)

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            val = planeRestData['Planes'][r][c]
            try:
                new = float(planeTable.GetValue(r,c))
                if new <= 0.:
                    raise ValueError
                planeRestData['Planes'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdatePlaneRestr,planeRestData)                
            
        def OnChangeEsd(event):
            rows = Planes.GetSelectedRows()
            if not rows:
                return
            Planes.ClearSelection()
            val = planeList[rows[0]][3]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for plane',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    planeRestData['Planes'][r][3] = parm
            dlg.Destroy()
            UpdatePlaneRestr(planeRestData)                
                                            
        def OnDeleteRestraint(event):
            rows = Planes.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                planeList.remove(planeList[row])
            UpdatePlaneRestr(planeRestData)                
            
        PlaneRestr.DestroyChildren()
        dataDisplay = wx.Panel(PlaneRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(PlaneRestr,planeRestData),0,wx.ALIGN_CENTER_VERTICAL)

        planeList = planeRestData['Planes']
        if len(planeList):
            table = []
            rowLabels = []
            bad = []
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) atom','calc','obs','esd']
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
                        table.append([name[:-3],calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            else:                                
                colLabels = ['atom+SymOp','calc','obs','esd']
                for i,[indx,ops,obs,esd] in enumerate(planeList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestPlane(XYZ,Amat)
                        name = ''
                        for a,atom in enumerate(atoms):
                            name += atom+'+ ('+ops[a]+'),'
                            if (a+1)%3 == 0:
                                name += '\n'
                        table.append([name[:-1],calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del planeList[ibad]
            planeTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Planes = G2gd.GSGrid(PlaneRestr)
            Planes.SetTable(planeTable, True)
            Planes.AutoSizeColumns(False)
            Planes.AutoSizeRows(False)
            for r in range(len(planeList)):
                for c in range(3):
                    Planes.SetReadOnly(r,c,True)
                    Planes.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Planes.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            Planes.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2gd.wxID_RESTCHANGEESD)
            mainSizer.Add(Planes,0,)
        else:
            mainSizer.Add(wx.StaticText(PlaneRestr,-1,'No plane restraints for this phase'),0,)

        PlaneRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50       #make room for tab
        PlaneRestr.SetSize(Size)
        PlaneRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)
#        G2frame.dataFrame.setSizePosLeft(Size)
    
    def UpdateChiralRestr(chiralRestData):

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            val = chiralRestData['Volumes'][r][c]
            try:
                new = float(volumeTable.GetValue(r,c))
                if new <= 0.:
                    raise ValueError
                chiralRestData['Volumes'][r][c] = new
            except ValueError:
                pass            
            wx.CallAfter(UpdateChiralRestr,chiralRestData)                
            
        def OnDeleteRestraint(event):
            rows = Volumes.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                volumeList.remove(volumeList[row])
            UpdateChiralRestr(chiralRestData)                
            
        def OnChangeValue(event):
            rows = Volumes.GetSelectedRows()
            if not rows:
                return
            Volumes.ClearSelection()
            val = volumeList[rows[0]][2]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new value for chiral volume',val,[0.,360.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    chiralRestData['Volumes'][r][2] = parm
            dlg.Destroy()
            UpdateChiralRestr(chiralRestData)                

        def OnChangeEsd(event):
            rows = Volumes.GetSelectedRows()
            if not rows:
                return
            Volumes.ClearSelection()
            val = volumeList[rows[0]][3]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for chiral volume',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    chiralRestData['Volumes'][r][3] = parm
            dlg.Destroy()
            UpdateChiralRestr(chiralRestData)                
                                            
        ChiralRestr.DestroyChildren()
        dataDisplay = wx.Panel(ChiralRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(ChiralRestr,chiralRestData),0,wx.ALIGN_CENTER_VERTICAL)

        volumeList = chiralRestData['Volumes']
        if len(volumeList):
            table = []
            rowLabels = []
            bad = []
            Types = [wg.GRID_VALUE_STRING,]+3*[wg.GRID_VALUE_FLOAT+':10,2',]
            if 'macro' in General['Type']:
                colLabels = ['(res) O (res) A (res) B (res) C','calc','obs','esd']
                for i,[indx,ops,obs,esd] in enumerate(volumeList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,0,4)
                        name = ''
                        for atom in atoms:
                            name += '('+atom[1]+atom[0].strip()+atom[2]+') '+atom[3]+' '
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        calc = G2mth.getRestChiral(XYZ,Amat)
                        table.append([name,calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            else:
                colLabels = ['O+SymOp  A+SymOp  B+SymOp  C+SymOp)','calc','obs','esd']
                for i,[indx,ops,obs,esd] in enumerate(volumeList):
                    try:
                        atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                        name = atoms[0]+'+('+ops[0]+') '+atoms[1]+'+('+ops[1]+') '+atoms[2]+ \
                            '+('+ops[2]+') '+atoms[3]+'+('+ops[3]+')'
                        XYZ = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cx,3))
                        XYZ = G2mth.getSyXYZ(XYZ,ops,SGData)
                        calc = G2mth.getRestChiral(XYZ,Amat)
                        table.append([name,calc,obs,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del volumeList[ibad]
            volumeTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Volumes = G2gd.GSGrid(ChiralRestr)
            Volumes.SetTable(volumeTable, True)
            Volumes.AutoSizeColumns(False)
            for r in range(len(volumeList)):
                for c in range(2):
                    Volumes.SetReadOnly(r,c,True)
                    Volumes.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Volumes.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            Volumes.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeValue, id=G2gd.wxID_RESRCHANGEVAL)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2gd.wxID_RESTCHANGEESD)
            mainSizer.Add(Volumes,0,)
        else:
            mainSizer.Add(wx.StaticText(ChiralRestr,-1,'No chiral volume restraints for this phase'),0,)

        ChiralRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50       #make room for tab
        ChiralRestr.SetSize(Size)
        ChiralRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)
    
    def UpdateTorsionRestr(torsionRestData):

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            val = torsionRestData['Torsions'][r][c]
            try:
                new = float(torsionTable.GetValue(r,c))
                if new <= 0. or new > 5.:
                    raise ValueError
                torsionRestData['Torsions'][r][3] = new     #only esd is editable
            except ValueError:
                pass            
            wx.CallAfter(UpdateTorsionRestr,torsionRestData)                
            
        def OnDeleteRestraint(event):
            rows = Torsions.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                torsionList.remove(torsionList[row])
            UpdateTorsionRestr(torsionRestData)                
            
        def OnChangeEsd(event):
            rows = Torsions.GetSelectedRows()
            if not rows:
                return
            Torsions.ClearSelection()
            val = torsionList[rows[0]][4]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for torsion restraints',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    torsionRestData['Torsions'][r][4] = parm
            dlg.Destroy()
            UpdateTorsionRestr(torsionRestData)                
                                            
        TorsionRestr.DestroyChildren()
        dataDisplay = wx.Panel(TorsionRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(TorsionRestr,torsionRestData),0,wx.ALIGN_CENTER_VERTICAL)
        
        coeffDict = torsionRestData['Coeff']
        torsionList = torsionRestData['Torsions']
        mainSizer.Add(wx.StaticText(TorsionRestr,-1,'Torsion restraints:'),0,wx.ALIGN_CENTER_VERTICAL)
        if len(torsionList):
            table = []
            rowLabels = []
            bad = []
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
                        table.append([name,cofName,tor,calc,restr,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del torsionList[ibad]
            torsionTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Torsions = G2gd.GSGrid(TorsionRestr)
            Torsions.SetTable(torsionTable, True)
            Torsions.AutoSizeColumns(False)
            for r in range(len(torsionList)):
                for c in range(5):
                    Torsions.SetReadOnly(r,c,True)
                    Torsions.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Torsions.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            Torsions.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2gd.wxID_RESTCHANGEESD)
            mainSizer.Add(Torsions,0,)
            
            mainSizer.Add((5,5))
            mainSizer.Add(wx.StaticText(TorsionRestr,-1,'Torsion function coefficients:'),0,wx.ALIGN_CENTER_VERTICAL)
            table = []
            rowLabels = []
            Types = 9*[wg.GRID_VALUE_FLOAT+':10,4',]
            colLabels = ['Mag A','Pos A','Width A','Mag B','Pos B','Width B','Mag C','Pos C','Width C']
            for item in coeffDict:
                rowLabels.append(item)
                table.append(coeffDict[item])
            coeffTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Coeff = G2gd.GSGrid(TorsionRestr)
            Coeff.SetTable(coeffTable, True)
            Coeff.AutoSizeColumns(False)
            for r in range(len(coeffDict)):
                for c in range(9):
                    Coeff.SetReadOnly(r,c,True)
                    Coeff.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            mainSizer.Add(Coeff,0,)
        else:
            mainSizer.Add(wx.StaticText(TorsionRestr,-1,'No torsion restraints for this phase'),0,)

        TorsionRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50       #make room for tab
        TorsionRestr.SetSize(Size)
        TorsionRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)

    def UpdateRamaRestr(ramaRestData):

        def OnCellChange(event):
            r,c =  event.GetRow(),event.GetCol()
            val = ramaRestData['Ramas'][r][c]
            try:
                new = float(ramaTable.GetValue(r,c))
                if new <= 0. or new > 5.:
                    raise ValueError
                ramaRestData['Ramas'][r][4] = new     #only esd is editable
            except ValueError:
                pass            
            wx.CallAfter(UpdateRamaRestr,ramaRestData)                
            
        def OnDeleteRestraint(event):
            rows = Ramas.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                ramaList.remove(ramaList[row])
            UpdateRamaRestr(ramaRestData)                
            
        def OnChangeEsd(event):
            rows = Ramas.GetSelectedRows()
            if not rows:
                return
            Ramas.ClearSelection()
            val = ramaList[rows[0]][4]
            dlg = G2phG.SingleFloatDialog(G2frame,'New value','Enter new esd for energy',val,[0.,5.],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    ramaRestData['Ramas'][r][4] = parm
            dlg.Destroy()
            UpdateRamaRestr(ramaRestData)                
                                            
        RamaRestr.DestroyChildren()
        dataDisplay = wx.Panel(RamaRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(RamaRestr,ramaRestData),0,wx.ALIGN_CENTER_VERTICAL)

        ramaList = ramaRestData['Ramas']
        coeffDict = ramaRestData['Coeff']
        if len(ramaList):
            mainSizer.Add(wx.StaticText(RamaRestr,-1,'Ramachandran restraints:'),0,wx.ALIGN_CENTER_VERTICAL)
            table = []
            rowLabels = []
            bad = []
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
                        table.append([name,cofName,phi,psi,calc,restr,esd])
                        rowLabels.append(str(i))
                    except KeyError:
                        print '**** WARNING - missing atom - restraint deleted ****'
                        bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del ramaList[ibad]
            ramaTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Ramas = G2gd.GSGrid(RamaRestr)
            Ramas.SetTable(ramaTable, True)
            Ramas.AutoSizeColumns(False)
            for r in range(len(ramaList)):
                for c in range(6):
                    Ramas.SetReadOnly(r,c,True)
                    Ramas.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            Ramas.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            Ramas.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeEsd, id=G2gd.wxID_RESTCHANGEESD)
            mainSizer.Add(Ramas,0,)
        else:
            mainSizer.Add(wx.StaticText(RamaRestr,-1,'No Ramachandran restraints for this phase'),0,)
        mainSizer.Add((5,5))
        mainSizer.Add(wx.StaticText(RamaRestr,-1,'Ramachandran function coefficients:'),0,wx.ALIGN_CENTER_VERTICAL)
        if len(coeffDict):
            table = []
            rowLabels = []
            Types = 6*[wg.GRID_VALUE_FLOAT+':10,4',]
            colLabels = ['Mag','Pos phi','Pos psi','sig(phi)','sig(psi)','sig(cov)']
            for item in coeffDict:
                for i,term in enumerate(coeffDict[item]):
                    rowLabels.append(item+' term:'+str(i))
                    table.append(term)
            coeffTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Coeff = G2gd.GSGrid(RamaRestr)
            Coeff.SetTable(coeffTable, True)
            Coeff.AutoSizeColumns(False)
            for r in range(Coeff.GetNumberRows()):
                for c in range(6):
                    Coeff.SetReadOnly(r,c,True)
                    Coeff.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
            mainSizer.Add(Coeff,0,)

        RamaRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50       #make room for tab
        RamaRestr.SetSize(Size)
        RamaRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)

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
            r = ChemComps.GetSelectedRows()[0]
            if not rows:
                return
            rowLabl = ChemComps.GetRowLabelValue(r)
            row = int(rowLabl.split(':')[1])
            if 'Restr' in rowLabl:
                del chemcompList[row]
            else:
                term = int(rowLabl.split(':')[2])
                del chemcompList[row][0][term]
                del chemcompList[row][1][term]
            UpdateChemcompRestr(chemcompRestData)                
            
        def OnChangeValue(event):
            rows = ChemComps.GetSelectedRows()
            if not rows:
                return
            ChemComps.ClearSelection()
            dlg = G2phG.SingleFloatDialog(G2frame,'New value',
                'Enter new value for restraint multiplier',1.0,[-1.e6,1.e6],'%.2f')
            if dlg.ShowModal() == wx.ID_OK:
                parm = dlg.GetValue()
                for r in rows:
                    rowLabl = ChemComps.GetRowLabelValue(r)
                    if 'term' in rowLabl:
                        items = rowLabl.split(':')
                        chemcompRestData['Sites'][int(items[1])][1][int(items[2])] = parm
            dlg.Destroy()
            UpdateChemcompRestr(chemcompRestData)                

        ChemCompRestr.DestroyChildren()
        dataDisplay = wx.Panel(ChemCompRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(ChemCompRestr,chemcompRestData),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(wx.StaticText(ChemCompRestr,-1, 
            'NB: The chemical restraint sum is over the unit cell contents'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)

        chemcompList = chemcompRestData['Sites']
        if len(chemcompList):
            table = []
            rowLabels = []
            bad = []
            Types = [wg.GRID_VALUE_STRING,]+5*[wg.GRID_VALUE_FLOAT+':10,2',]
            colLabels = ['Atoms','mul*frac','factor','calc','obs','esd']
            for i,[indx,factors,obs,esd] in enumerate(chemcompList):
                try:
                    atoms = G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,ct-1)
                    mul = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs+1))
                    frac = np.array(G2mth.GetAtomItemsById(Atoms,AtLookUp,indx,cs-1))
                    mulfrac = mul*frac
                    calcs = mul*frac*factors
                    for iatm,[atom,mf,fr,clc] in enumerate(zip(atoms,mulfrac,factors,calcs)):
                        table.append([atom,mf,fr,clc,'',''])
                        rowLabels.append('term:'+str(i)+':'+str(iatm))
                    table.append(['Sum','','',np.sum(calcs),obs,esd])
                    rowLabels.append('Restr:'+str(i))
                except KeyError:
                    print '**** WARNING - missing atom - restraint deleted ****'
                    bad.append(i)
            if len(bad):
                bad.reverse()
                for ibad in bad:
                    del chemcompList[ibad]
            chemcompTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            ChemComps = G2gd.GSGrid(ChemCompRestr)
            ChemComps.SetTable(chemcompTable, True)
            ChemComps.AutoSizeColumns(False)
            for r in range(chemcompTable.GetNumberRows()):
                for c in range(2):
                    ChemComps.SetReadOnly(r,c,True)
                    ChemComps.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                if 'Restr' in ChemComps.GetRowLabelValue(r):
                    for c in range(4):
                        ChemComps.SetReadOnly(r,c,True)
                        ChemComps.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                    for c in [1,2]:
                        ChemComps.SetCellTextColour(r,c,VERY_LIGHT_GREY)
                else:
                    for c in [3,4,5]:
                        ChemComps.SetReadOnly(r,c,True)
                        ChemComps.SetCellStyle(r,c,VERY_LIGHT_GREY,True)
                    for c in [4,5]:
                        ChemComps.SetCellTextColour(r,c,VERY_LIGHT_GREY)
                        
            ChemComps.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK,OnRowSelect)
            ChemComps.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnChangeValue, id=G2gd.wxID_RESRCHANGEVAL)
            mainSizer.Add(ChemComps,0,)
        else:
            mainSizer.Add(wx.StaticText(ChemCompRestr,-1,'No chemical composition restraints for this phase'),0,)

        ChemCompRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50       #make room for tab
        ChemCompRestr.SetSize(Size)
        ChemCompRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)
    
        
    def UpdateTextureRestr(textureRestData):
            
        def OnDeleteRestraint(event):
            rows = Textures.GetSelectedRows()
            if not rows:
                return
            rows.sort()
            rows.reverse()
            for row in rows:
                textureList.remove(textureList[row])
            wx.CallAfter(UpdateTextureRestr,textureRestData)                
            
        def OnCellChange(event):
            r,c = event.GetRow(),event.GetCol()
            val = textureRestData['HKLs'][r][c]
            try:
                if c == 1:  #grid size
                    new = int(textureTable.GetValue(r,c))
                    if new < 6 or new > 24:
                        raise valueError
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

        TextureRestr.DestroyChildren()
        dataDisplay = wx.Panel(TextureRestr)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add((5,5),0)
        mainSizer.Add(WtBox(TextureRestr,textureRestData),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(wx.StaticText(TextureRestr,-1, 
            'NB: The texture restraints suppress negative pole figure values for the selected HKLs\n'
            '    "unit esd" gives a bias toward a flatter polefigure'),0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add((5,5),0)

        textureList = textureRestData['HKLs']
        textureData = General['SH Texture']
        if len(textureList):
            table = []
            rowLabels = []
            Types = [wg.GRID_VALUE_STRING,wg.GRID_VALUE_LONG,wg.GRID_VALUE_FLOAT+':10,2',
                wg.GRID_VALUE_BOOL,wg.GRID_VALUE_FLOAT+':10,2']
            colLabels = ['HKL','grid','neg esd','use unit?','unit esd']
            for i,[hkl,grid,esd1,ifesd2,esd2] in enumerate(textureList):
                table.append(['%d %d %d'%(hkl[0],hkl[1],hkl[2]),grid,esd1,ifesd2,esd2])
                rowLabels.append(str(i))
            textureTable = G2gd.Table(table,rowLabels=rowLabels,colLabels=colLabels,types=Types)
            Textures = G2gd.GSGrid(TextureRestr)
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
            Textures.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
            G2frame.dataFrame.Bind(wx.EVT_MENU, OnDeleteRestraint, id=G2gd.wxID_RESTDELETE)
            mainSizer.Add(Textures,0,)
        else:
            mainSizer.Add(wx.StaticText(TextureRestr,-1,'No texture restraints for this phase'),0,)
        TextureRestr.SetSizer(mainSizer)
        Size = mainSizer.Fit(G2frame.dataFrame)
        Size[0] = 600
        Size[1] += 50       #make room for tab
        TextureRestr.SetSize(Size)
        TextureRestr.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        G2frame.dataFrame.SetSize(Size)
            
    def OnPageChanged(event):
        page = event.GetSelection()
        text = G2frame.dataDisplay.GetPageText(page)
        G2frame.dataFrame.RestraintEdit.SetLabel(G2gd.wxID_RESRCHANGEVAL,'Change value')
        SetStatusLine('')
        if text == 'Bond restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,True)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,True)
            bondRestData = restrData['Bond']
            UpdateBondRestr(bondRestData)
        elif text == 'Angle restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,True)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,True)
            angleRestData = restrData['Angle']
            UpdateAngleRestr(angleRestData)
        elif text == 'Plane restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,True)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,False)
            planeRestData = restrData['Plane']
            UpdatePlaneRestr(planeRestData)
        elif text == 'Chiral restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,False)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,True)
            chiralRestData = restrData['Chiral']
            UpdateChiralRestr(chiralRestData)
        elif text == 'Torsion restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,False)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,False)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_AARESTRAINTPLOT,True)
            torsionRestData = restrData['Torsion']
            UpdateTorsionRestr(torsionRestData)
        elif text == 'Ramachandran restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,False)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,False)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_AARESTRAINTPLOT,True)
            ramaRestData = restrData['Rama']
            UpdateRamaRestr(ramaRestData)
            G2plt.PlotRama(G2frame,phaseName,rama,ramaName)
        elif text == 'Chem. comp. restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,True)
            G2frame.dataFrame.RestraintEdit.SetLabel(G2gd.wxID_RESRCHANGEVAL,'Change factor')
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,True)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTCHANGEESD,False)
            chemcompRestData = restrData['ChemComp']
            UpdateChemcompRestr(chemcompRestData)
        elif text == 'Texture restraints':
            G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTRAINTADD,True)
            G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESRCHANGEVAL,True)
            textureRestData = restrData['Texture']
            UpdateTextureRestr(textureRestData)
            
        event.Skip()

    def SetStatusLine(text):
        Status.SetStatusText(text)                                      
        
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
        
    G2gd.SetDataMenuBar(G2frame,G2frame.dataFrame.RestraintMenu)
    G2frame.dataFrame.SetLabel('restraints for '+phaseName)
    if not G2frame.dataFrame.GetStatusBar():
        Status = G2frame.dataFrame.CreateStatusBar()
    SetStatusLine('')
    
    G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTSELPHASE,False)
    if len(Phases) > 1:
        G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_RESTSELPHASE,True)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnSelectPhase, id=G2gd.wxID_RESTSELPHASE)
    G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddRestraint, id=G2gd.wxID_RESTRAINTADD)
    if 'macro' in phasedata['General']['Type']:
        G2frame.dataFrame.RestraintEdit.Enable(G2gd.wxID_AARESTRAINTADD,True)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnAddAARestraint, id=G2gd.wxID_AARESTRAINTADD)
        G2frame.dataFrame.Bind(wx.EVT_MENU, OnPlotAARestraint, id=G2gd.wxID_AARESTRAINTPLOT)
    G2frame.dataDisplay = G2gd.GSNoteBook(parent=G2frame.dataFrame,size=G2frame.dataFrame.GetClientSize())
    
    BondRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(BondRestr,'Bond restraints')
    AngleRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(AngleRestr,'Angle restraints')
    PlaneRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(PlaneRestr,'Plane restraints')
    ChiralRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(ChiralRestr,'Chiral restraints')
    if 'macro' in General['Type']:
        TorsionRestr = wx.ScrolledWindow(G2frame.dataDisplay)
        G2frame.dataDisplay.AddPage(TorsionRestr,'Torsion restraints')
        RamaRestr = wx.ScrolledWindow(G2frame.dataDisplay)
        G2frame.dataDisplay.AddPage(RamaRestr,'Ramachandran restraints')
    ChemCompRestr = wx.ScrolledWindow(G2frame.dataDisplay)
    G2frame.dataDisplay.AddPage(ChemCompRestr,'Chem. comp. restraints')
    if General['SH Texture']['Order']:
        TextureRestr = wx.ScrolledWindow(G2frame.dataDisplay)
        G2frame.dataDisplay.AddPage(TextureRestr,'Texture restraints')
    
    UpdateBondRestr(restrData['Bond'])

    G2frame.dataDisplay.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, OnPageChanged)
