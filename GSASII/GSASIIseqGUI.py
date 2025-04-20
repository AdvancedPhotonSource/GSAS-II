# -*- coding: utf-8 -*-
#GSASIIseqGUI - Sequential Results Display routines
'''
Routines for Sequential Results & Cluster Analysis dataframes follow.
'''
from __future__ import division, print_function
import platform
import copy
import re
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import scipy.optimize as so
try:
    import wx
    import wx.grid as wg
except ImportError:
    pass
from . import GSASIIpath
from . import GSASIImath as G2mth
from . import GSASIImiscGUI as G2IO
from . import GSASIIdataGUI as G2gd
from . import GSASIIstrIO as G2stIO
from . import GSASIIlattice as G2lat
from . import GSASIIplot as G2plt
from . import GSASIIpwdplot as G2pwpl
from . import GSASIImapvars as G2mv
from . import GSASIIobj as G2obj
from . import GSASIIexprGUI as G2exG
from . import GSASIIctrlGUI as G2G
WACV = wx.ALIGN_CENTER_VERTICAL

#####  Display of Sequential Results ##########################################
def UpdateSeqResults(G2frame,data,prevSize=None):
    """
    Called when any data tree entry is selected that has 'Sequential' in the name
    to show results from any sequential analysis.

    :param wx.Frame G2frame: main GSAS-II data tree windows

    :param dict data: a dictionary containing the following items:

            * 'histNames' - list of histogram names in order as processed by Sequential Refinement
            * 'varyList' - list of variables - identical over all refinements in sequence
              note that this is the original list of variables, prior to processing
              constraints.
            * 'variableLabels' -- a dict of labels to be applied to each parameter
              (this is created as an empty dict if not present in data).
            * keyed by histName - dictionaries for all data sets processed, which contains:

              * 'variables'- result[0] from leastsq call
              * 'varyList' - list of variables passed to leastsq call (not same as above)
              * 'sig' - esds for variables
              * 'covMatrix' - covariance matrix from individual refinement
              * 'title' - histogram name; same as dict item name
              * 'newAtomDict' - new atom parameters after shifts applied
              * 'newCellDict' - refined cell parameters after shifts to A0-A5 from Dij terms applied'
    """
    def GetSampleParms():
        '''Make a dictionary of the sample parameters that are not the same over the
        refinement series. Controls here is local
        '''
        if 'IMG' in histNames[0]:
            sampleParmDict = {'Sample load':[],}
        elif 'PDF' in histNames[0]:
            sampleParmDict = {'Temperature':[]}
        else:
            sampleParmDict = {'Temperature':[],'Pressure':[],'Time':[],
                'FreePrm1':[],'FreePrm2':[],'FreePrm3':[],'Omega':[],
                'Chi':[],'Phi':[],'Azimuth':[],}
        Controls = G2frame.GPXtree.GetItemPyData(
            G2gd.GetGPXtreeItemId(G2frame,G2frame.root, 'Controls'))
        sampleParm = {}
        for name in histNames:
            if 'IMG' in name or 'PDF' in name:
                if name not in data:
                    continue
                for item in sampleParmDict:
                    sampleParmDict[item].append(data[name]['parmDict'].get(item,0))
            else:
                if 'PDF' in name:
                    name = 'PWDR' + name[4:]
                Id = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                if Id:
                    sampleData = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,Id,'Sample Parameters'))
                else:  # item missing from tree! stick in NaN's!
                    sampleData = {}
                for item in sampleParmDict:
                    sampleParmDict[item].append(sampleData.get(item,np.nan))
        for item in sampleParmDict:
            if sampleParmDict[item]:
                frstValue = sampleParmDict[item][0]
                if np.any(np.array(sampleParmDict[item])-frstValue):
                    if item.startswith('FreePrm'):
                        sampleParm[Controls[item]] = sampleParmDict[item]
                    else:
                        sampleParm[item] = sampleParmDict[item]
        return sampleParm

    def GetColumnInfo(col):
        '''returns column label, lists of values and errors (or None) for each column in the table
        for plotting. The column label is reformatted from Unicode to MatPlotLib encoding
        '''
        colName = G2frame.SeqTable.GetColLabelValue(col)
        plotName = variableLabels.get(colName,colName)
        plotName = plotSpCharFix(plotName)
        return plotName,G2frame.colList[col],G2frame.colSigs[col]

    def PlotSelectedColRow(calltyp='',event=None):
        '''Called to plot a selected column or row by clicking
        on a row or column label. N.B. This is called for LB click
        after the event is processed so that the column or row has been
        selected. For RB clicks, the event is processed here

        :param str calltyp: ='left'/'right', specifies if this was
          a left- or right-click, where a left click on row
          plots histogram; Right click on row plots V-C matrix;
          Left or right click on column: plots values in column
        :param obj event: from  wx.EVT_GRID_LABEL_RIGHT_CLICK
        '''
        rows = []
        cols = []
        if calltyp == 'left':
            cols = G2frame.dataDisplay.GetSelectedCols()
            rows = G2frame.dataDisplay.GetSelectedRows()
        elif calltyp == 'right':
            r,c = event.GetRow(),event.GetCol()
            if r > -1:
                rows = [r,]
            elif c > -1:
                cols =  [c,]
        if cols and calltyp == 'left':
            G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
        elif cols and calltyp == 'right':
            SetLabelString(cols[0])  #only the 1st one selected!
        elif rows and calltyp == 'left':
            name = histNames[rows[0]]       #only does 1st one selected
            if name.startswith('PWDR'):
                pickId = G2frame.PickId
                G2frame.PickId = G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, name)
                G2pwpl.PlotPatterns(G2frame,newPlot=False,plotType='PWDR')
                G2frame.PickId = pickId
            elif name.startswith('PDF'):
                pickId = G2frame.PickId
                G2frame.PickId = G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame, G2frame.root, name)
                PFdata = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.PickId,'PDF Controls'))
                G2plt.PlotISFG(G2frame,PFdata,plotType='G(R)')
                G2frame.PickId = pickId
            else:
                return
        elif rows and calltyp == 'right':
            name = histNames[rows[0]]       #only does 1st one selected
            if name.startswith('PWDR'):
                if len(data[name].get('covMatrix',[])):
                    G2plt.PlotCovariance(G2frame,data[name])
            elif name.startswith('PDF'):
                print('make structure plot')
        else:
            G2frame.ErrorDialog(
                'Select row or columns',
                'Nothing selected in table. Click on column or row label(s) to plot. N.B. Grid selection can be a bit funky.'
                )

    def SetLabelString(col):
        '''Define or edit the label for a column in the table, to be used
        as a tooltip and for plotting
        '''
        if col < 0 or col > len(colLabels):
            return
        var = colLabels[col]
        lbl = variableLabels.get(var,G2obj.fmtVarDescr(var))
        head = u'Set a new name for variable {} (column {})'.format(var,col)
        dlg = G2G.SingleStringDialog(G2frame,'Set variable label',
                                 head,lbl,size=(400,-1))
        if dlg.Show():
            variableLabels[var] = dlg.GetValue()
            dlg.Destroy()
            wx.CallAfter(UpdateSeqResults,G2frame,data) # redisplay variables
        else:
            dlg.Destroy()

    def PlotLeftSelect(event):
        'Called by a left MB click on a row or column label. '
        event.Skip()
        wx.CallAfter(PlotSelectedColRow,'left')

    def PlotRightSelect(event):
        'Called by a right MB click on a row or column label'
        PlotSelectedColRow('right',event)

    def OnPlotSelSeq(event):
        'plot the selected columns or row from menu command'
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        rows = G2frame.dataDisplay.GetSelectedRows()
        if cols:
            G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)
        elif rows:
            name = histNames[rows[0]]       #only does 1st one selected
            G2plt.PlotCovariance(G2frame,data[name])
        else:
            G2frame.ErrorDialog('Select columns',
                'No columns or rows selected in table. Click on row or column labels to select fields for plotting.')

    def OnAveSelSeq(event):
        'average the selected columns from menu command'
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        useCol =  ~np.array(G2frame.SeqTable.GetColValues(1),dtype=bool)
        if cols:
            for col in cols:
                items = GetColumnInfo(col)[1]
                noneMask = np.array([item is None for item in items])
                info = ma.array(items,mask=useCol+noneMask)
                ave = ma.mean(ma.compressed(info))
                sig = ma.std(ma.compressed(info))
                print (u' Average for '+G2frame.SeqTable.GetColLabelValue(col)+u': '+'%.6g'%(ave)+u' +/- '+u'%.6g'%(sig))
        else:
            G2frame.ErrorDialog('Select columns',
                'No columns selected in table. Click on column labels to select fields for averaging.')

    def OnSelectUse(event):
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select rows to use','Select rows',histNames)
        sellist = [i for i,item in enumerate(G2frame.colList[1]) if item]
        dlg.SetSelections(sellist)
        if dlg.ShowModal() == wx.ID_OK:
            sellist = dlg.GetSelections()
            for row in range(G2frame.SeqTable.GetNumberRows()):
                G2frame.SeqTable.SetValue(row,1,False)
                G2frame.colList[1][row] = False
            for row in sellist:
                G2frame.SeqTable.SetValue(row,1,True)
                G2frame.colList[1][row] = True
            G2frame.dataDisplay.ForceRefresh()
        dlg.Destroy()

    def OnRenameSelSeq(event):
        cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        colNames = [G2frame.SeqTable.GetColLabelValue(c) for c in cols]
        newNames = colNames[:]
        for i,name in enumerate(colNames):
            if name in variableLabels:
                newNames[i] = variableLabels[name]
        if not cols:
            G2frame.ErrorDialog('Select columns',
                'No columns selected in table. Click on column labels to select fields for rename.')
            return
        dlg = G2G.MultiStringDialog(G2frame.dataDisplay,'Set column names',colNames,newNames)
        if dlg.Show():
            newNames = dlg.GetValues()
            variableLabels.update(dict(zip(colNames,newNames)))
        data['variableLabels'] = variableLabels
        dlg.Destroy()
        UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
        G2plt.PlotSelectedSequence(G2frame,cols,GetColumnInfo,SelectXaxis)

    def OnSaveSelSeqCSV(event):
        'export the selected columns to a .csv file from menu command'
        OnSaveSelSeq(event,csv=True)

    def OnSaveSeqCSV(event):
        'export all columns to a .csv file from menu command'
        OnSaveSelSeq(event,csv=True,allcols=True)

    def OnSaveSelSeq(event,csv=False,allcols=False):
        'export the selected columns to a .txt or .csv file from menu command'
        def WriteLine(line):
            if '2' in platform.python_version_tuple()[0]:
                SeqFile.write(G2obj.StripUnicode(line))
            else:
                SeqFile.write(line)

        def WriteCSV():
            def WriteList(headerItems):
                line = ''
                for lbl in headerItems:
                    if line: line += ','
                    line += '"'+lbl+'"'
                return line
            head = ['name']
            for col in cols:
                # Excel does not like unicode
                item = G2obj.StripUnicode(G2frame.SeqTable.GetColLabelValue(col))
                if col in havesig:
                    head += [item,'esd-'+item]
                else:
                    head += [item]
            WriteLine(WriteList(head)+'\n')
            for row,name in enumerate(saveNames):
                line = '"'+saveNames[row]+'"'
                for col in cols:
                    if saveData[col][row] is None:
                        if col in havesig:
#                            line += ',0.0,0.0'
                            line += ',,'
                        else:
#                            line += ',0.0'
                            line += ','
                    else:
                        if col in havesig:
                            line += ','+str(saveData[col][row])+','+str(saveSigs[col][row])
                        else:
                            line += ','+str(saveData[col][row])
                WriteLine(line+'\n')
        def WriteSeq():
            lenName = len(saveNames[0])
            line = '  %s  '%('name'.center(lenName))
            for col in cols:
                item = G2frame.SeqTable.GetColLabelValue(col)
                if col in havesig:
                    line += ' %12s %12s '%(item.center(12),'esd'.center(12))
                else:
                    line += ' %12s '%(item.center(12))
            WriteLine(line+'\n')
            for row,name in enumerate(saveNames):
                line = " '%s' "%(saveNames[row])
                for col in cols:
                    if col in havesig:
                        try:
                            line += ' %12.6f %12.6f '%(saveData[col][row],saveSigs[col][row])
                        except TypeError:
                            line += '                           '
                    else:
                        try:
                            line += ' %12.6f '%saveData[col][row]
                        except TypeError:
                            line += '              '
                WriteLine(line+'\n')

        # start of OnSaveSelSeq code
        if allcols:
            cols = range(G2frame.SeqTable.GetNumberCols())
        else:
            cols = sorted(G2frame.dataDisplay.GetSelectedCols()) # ignore selection order
        nrows = G2frame.SeqTable.GetNumberRows()
        if not cols:
            choices = [G2frame.SeqTable.GetColLabelValue(r) for r in range(G2frame.SeqTable.GetNumberCols())]
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select columns to write',
                'select columns',choices)
            #dlg.SetSelections()
            if dlg.ShowModal() == wx.ID_OK:
                cols = dlg.GetSelections()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
            #G2frame.ErrorDialog('Select columns',
            #                 'No columns selected in table. Click on column labels to select fields for output.')
            #return
        saveNames = [G2frame.SeqTable.GetRowLabelValue(r) for r in range(nrows)]
        saveData = {}
        saveSigs = {}
        havesig = []
        for col in cols:
            name,vals,sigs = GetColumnInfo(col)
            saveData[col] = vals
            if sigs:
                havesig.append(col)
                saveSigs[col] = sigs
        if csv:
            wild = 'CSV output file (*.csv)|*.csv'
        else:
            wild = 'Text output file (*.txt)|*.txt'
        pth = G2G.GetExportPath(G2frame)
        dlg = wx.FileDialog(
            G2frame,
            'Choose text output file for your selection', pth, '',
            wild,wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                SeqTextFile = dlg.GetPath()
                SeqTextFile = G2IO.FileDlgFixExt(dlg,SeqTextFile)
                SeqFile = open(SeqTextFile,'w')
                if csv:
                    WriteCSV()
                else:
                    WriteSeq()
                SeqFile.close()
        finally:
            dlg.Destroy()

    def striphist(var,insChar=''):
        'strip a histogram number from a var name'
        sv = var.split(':')
        if len(sv) <= 1: return var
        if sv[1]:
            sv[1] = insChar
        return ':'.join(sv)

    def plotSpCharFix(lbl):
        'Change selected unicode characters to their matplotlib equivalent'
        for u,p in [
            (u'\u03B1',r'$\alpha$'),
            (u'\u03B2',r'$\beta$'),
            (u'\u03B3',r'$\gamma$'),
            (u'\u0394\u03C7',r'$\Delta\chi$'),
            ]:
            lbl = lbl.replace(u,p)
        return lbl

    def SelectXaxis():
        'returns a selected column number (or None) as the X-axis selection'
        ncols = G2frame.SeqTable.GetNumberCols()
        colNames = [G2frame.SeqTable.GetColLabelValue(r) for r in range(ncols)]
        dlg = G2G.G2SingleChoiceDialog(
            G2frame.dataDisplay,
            'Select x-axis parameter for\nplot (Cancel=sequence #)',
            'Select X-axis',
            colNames)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                col = dlg.GetSelection()
            else:
                col = None
        finally:
            dlg.Destroy()
        return col

    def EnablePseudoVarMenus():
        'Enables or disables the PseudoVar menu items that require existing defs'
        if data['SeqPseudoVars']:
            val = True
        else:
            val = False
        G2frame.dataWindow.SequentialPvars.Enable(G2G.wxID_DELSEQVAR,val)
        G2frame.dataWindow.SequentialPvars.Enable(G2G.wxID_EDITSEQVAR,val)

    def DelPseudoVar(event):
        'Ask the user to select a pseudo var expression to delete'
        choices = list(data['SeqPseudoVars'].keys())
        selected = G2G.ItemSelector(
            choices,G2frame,
            multiple=True,
            title='Select expressions to remove',
            header='Delete expression')
        if selected is None: return
        for item in selected:
            del data['SeqPseudoVars'][choices[item]]
        if selected:
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def EditPseudoVar(event):
        'Edit an existing pseudo var expression'
        choices = list(data['SeqPseudoVars'].keys())
        if len(choices) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                choices,G2frame,
                multiple=False,
                title='Select an expression to edit',
                header='Edit expression')
        if selected is not None:
            dlg = G2exG.ExpressionDialog(
                G2frame.dataDisplay,PSvarDict,
                data['SeqPseudoVars'][choices[selected]],
                header="Edit the PseudoVar expression",
                VarLabel="PseudoVar #"+str(selected+1),
                fit=False)
            newobj = dlg.Show(True)
            if newobj:
                calcobj = G2obj.ExpressionCalcObj(newobj)
                del data['SeqPseudoVars'][choices[selected]]
                data['SeqPseudoVars'][calcobj.eObj.expression] = newobj
                UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def AddNewPseudoVar(event):
        'Create a new pseudo var expression'
        dlg = G2exG.ExpressionDialog(G2frame.dataDisplay,PSvarDict,
            header='Enter an expression for a PseudoVar here',
            VarLabel = "New PseudoVar",fit=False)
        obj = dlg.Show(True)
        dlg.Destroy()
        if obj:
            calcobj = G2obj.ExpressionCalcObj(obj)
            data['SeqPseudoVars'][calcobj.eObj.expression] = obj
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def AddNewDistPseudoVar(event):
        obj = None
        dlg = G2exG.BondDialog(
            G2frame.dataDisplay,Phases,PSvarDict,
            VarLabel = "New Bond")
        if dlg.ShowModal() == wx.ID_OK:
            pName,Oatom,Tatom = dlg.GetSelection()
            if Tatom:
                Phase = Phases[pName]
                General = Phase['General']
                cx,ct = General['AtomPtrs'][:2]
                pId = Phase['pId']
                SGData = General['SGData']
                sB = Tatom.find('(')+1
                symNo = 0
                if sB:
                    sF = Tatom.find(')')
                    symNo = int(Tatom[sB:sF])
                cellNo = [0,0,0]
                cB = Tatom.find('[')
                if cB>0:
                    cF = Tatom.find(']')+1
                    cellNo = eval(Tatom[cB:cF])
                Atoms = Phase['Atoms']
                aNames = [atom[ct-1] for atom in Atoms]
                oId = aNames.index(Oatom)
                tId = aNames.index(Tatom.split(' +')[0])
                # create an expression object
                obj = G2obj.ExpressionObj()
                obj.expression = 'Dist(%s,\n%s)'%(Oatom,Tatom.split(' d=')[0].replace(' ',''))
                obj.distance_dict = {'pId':pId,'SGData':SGData,'symNo':symNo,'cellNo':cellNo}
                obj.distance_atoms = [oId,tId]
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        if obj:
            data['SeqPseudoVars'][obj.expression] = obj
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def AddNewAnglePseudoVar(event):
        obj = None
        dlg = G2exG.AngleDialog(
            G2frame.dataDisplay,Phases,PSvarDict,
            header='Enter an Angle here',
            VarLabel = "New Angle")
        if dlg.ShowModal() == wx.ID_OK:
            pName,Oatom,Tatoms = dlg.GetSelection()
            if Tatoms:
                obj = G2obj.makeAngleObj(Phases[pName],Oatom,Tatoms)
        else:
            dlg.Destroy()
            return
        dlg.Destroy()
        if obj:
            data['SeqPseudoVars'][obj.expression] = obj
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables

    def UpdateParmDict(parmDict):
        '''generate the atom positions and the direct & reciprocal cell values,
        because they might be needed to evaluate the pseudovar
        #TODO - effect of ISO modes?
        '''
        Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],
                         ['A'+str(i) for i in range(6)])
                     )
        delList = []
        phaselist = []
        for item in parmDict:
            if ':' not in item: continue
            key = item.split(':')
            if len(key) < 3: continue
            # remove the dA[xyz] terms, they would only bring confusion
            if key[0] and key[0] not in phaselist: phaselist.append(key[0])
            if key[2].startswith('dA'):
                delList.append(item)
            # compute and update the corrected reciprocal cell terms using the Dij values
            elif key[2] in Ddict:
                akey = key[0]+'::'+Ddict[key[2]]
                parmDict[akey] += parmDict[item]
                delList.append(item)
        for item in delList:
            del parmDict[item]
        for i in phaselist:
            pId = int(i)
            # apply cell symmetry
            try:
                A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata[pId],parmDict,zeroDict[pId])
                # convert to direct cell & add the unique terms to the dictionary
                try:
                    dcell = G2lat.A2cell(A)
                except:
                    print('phase',pId,'Invalid cell tensor',A)
                    raise ValueError('Invalid cell tensor in phase '+str(pId))
                for i,val in enumerate(dcell):
                    if i in uniqCellIndx[pId]:
                        lbl = str(pId)+'::'+G2lat.cellUlbl[i]
                        parmDict[lbl] = val
                lbl = str(pId)+'::'+'Vol'
                parmDict[lbl] = G2lat.calc_V(A)
            except KeyError:
                pass
        return parmDict

    def EvalPSvarDeriv(calcobj,parmDict,sampleDict,var,ESD):
        '''Evaluate an expression derivative with respect to a
        GSAS-II variable name.

        Note this likely could be faster if the loop over calcobjs were done
        inside after the Dict was created.
        '''
        if not ESD:
            return 0.
        step = ESD/10
        Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],
                         ['A'+str(i) for i in range(6)])
                     )
        results = []
        phaselist = []
        VparmDict = sampleDict.copy()
        for incr in step,-step:
            VparmDict.update(parmDict.copy())
            # as saved, the parmDict has updated 'A[xyz]' values, but 'dA[xyz]'
            # values are not zeroed: fix that!
            VparmDict.update({item:0.0 for item in parmDict if 'dA' in item})
            VparmDict[var] += incr
            G2mv.Dict2Map(VparmDict) # apply constraints
            # generate the atom positions and the direct & reciprocal cell values now, because they might
            # needed to evaluate the pseudovar
            for item in VparmDict:
                if item in sampleDict:
                    continue
                if ':' not in item: continue
                key = item.split(':')
                if len(key) < 3: continue
                # apply any new shifts to atom positions
                if key[2].startswith('dA'):
                    VparmDict[''.join(item.split('d'))] += VparmDict[item]
                    VparmDict[item] = 0.0
                # compute and update the corrected reciprocal cell terms using the Dij values
                if key[2] in Ddict:
                    if key[0] not in phaselist: phaselist.append(key[0])
                    akey = key[0]+'::'+Ddict[key[2]]
                    VparmDict[akey] += VparmDict[item]
            for i in phaselist:
                pId = int(i)
                # apply cell symmetry
                try:
                    A,zeros = G2stIO.cellFill(str(pId)+'::',SGdata[pId],VparmDict,zeroDict[pId])
                    # convert to direct cell & add the unique terms to the dictionary
                    for i,val in enumerate(G2lat.A2cell(A)):
                        if i in uniqCellIndx[pId]:
                            lbl = str(pId)+'::'+G2lat.cellUlbl[i]
                            VparmDict[lbl] = val
                    lbl = str(pId)+'::'+'Vol'
                    VparmDict[lbl] = G2lat.calc_V(A)
                except KeyError:
                    pass
            # dict should be fully updated, use it & calculate
            calcobj.SetupCalc(VparmDict)
            results.append(calcobj.EvalExpression())
        if None in results:
            return None
        return (results[0] - results[1]) / (2.*step)

    def EnableParFitEqMenus():
        'Enables or disables the Parametric Fit menu items that require existing defs'
        if data['SeqParFitEqList']:
            val = True
        else:
            val = False
        G2frame.dataWindow.SequentialPfit.Enable(G2G.wxID_DELPARFIT,val)
        G2frame.dataWindow.SequentialPfit.Enable(G2G.wxID_EDITPARFIT,val)
        G2frame.dataWindow.SequentialPfit.Enable(G2G.wxID_DOPARFIT,val)

    def ParEqEval(Values,calcObjList,varyList):
        '''Evaluate the parametric expression(s)
        :param list Values: a list of values for each variable parameter
        :param list calcObjList: a list of :class:`GSASIIobj.ExpressionCalcObj`
          expression objects to evaluate
        :param list varyList: a list of variable names for each value in Values
        '''
        result = []
        for calcobj in calcObjList:
            calcobj.UpdateVars(varyList,Values)
            if calcobj.depSig:
                result.append((calcobj.depVal-calcobj.EvalExpression())/calcobj.depSig)
            else:
                result.append(calcobj.depVal-calcobj.EvalExpression())
        return result

    def DoParEqFit(event,eqObj=None):
        'Parametric fit minimizer'
        varyValueDict = {} # dict of variables and their initial values
        calcObjList = [] # expression objects, ready to go for each data point
        if eqObj is not None:
            eqObjList = [eqObj,]
        else:
            eqObjList = data['SeqParFitEqList']
        UseFlags = G2frame.SeqTable.GetColValues(1)
        for obj in eqObjList:
            # assemble refined vars for this equation
            varyValueDict.update({var:val for var,val in obj.GetVariedVarVal()})
            # lookup dependent var position
            depVar = obj.GetDepVar()
            if depVar in colLabels:
                indx = colLabels.index(depVar)
            else:
                raise Exception('Dependent variable '+depVar+' not found')
            # assemble a list of the independent variables
            indepVars = obj.GetIndependentVars()
            # loop over each datapoint
            for j,row in enumerate(zip(*G2frame.colList)):
                if not UseFlags[j]: continue
                # assemble equations to fit
                calcobj = G2obj.ExpressionCalcObj(obj)
                # prepare a dict of needed independent vars for this expression
                indepVarDict = {var:row[i] for i,var in enumerate(colLabels) if var in indepVars}
                calcobj.SetupCalc(indepVarDict)
                # values and sigs for current value of dependent var
                if row[indx] is None: continue
                calcobj.depVal = row[indx]
                if G2frame.colSigs[indx]:
                    calcobj.depSig = G2frame.colSigs[indx][j]
                else:
                    calcobj.depSig = 1.
                calcObjList.append(calcobj)
        # varied parameters
        varyList = varyValueDict.keys()
        values = varyValues = [varyValueDict[key] for key in varyList]
        if not varyList:
            print ('no variables to refine!')
            return
        try:
            result = so.leastsq(ParEqEval,varyValues,full_output=True,   #ftol=Ftol,
                args=(calcObjList,varyList))
            values = result[0]
            covar = result[1]
            if covar is None:
                raise Exception
            chisq = np.sum(result[2]['fvec']**2)
            GOF = np.sqrt(chisq/(len(calcObjList)-len(varyList)))
            esdDict = {}
            for i,avar in enumerate(varyList):
                esdDict[avar] = np.sqrt(covar[i,i])
        except:
            print('====> Fit failed')
            return
        print('==== Fit Results ====')
        print ('  chisq =  %.2f, GOF = %.2f'%(chisq,GOF))
        for obj in eqObjList:
            obj.UpdateVariedVars(varyList,values)
            ind = '      '
            print(u'  '+obj.GetDepVar()+u' = '+obj.expression)
            for var in obj.assgnVars:
                print(ind+var+u' = '+obj.assgnVars[var])
            for var in obj.freeVars:
                avar = "::"+obj.freeVars[var][0]
                val = obj.freeVars[var][1]
                if obj.freeVars[var][2]:
                    print(ind+var+u' = '+avar + " = " + G2mth.ValEsd(val,esdDict[avar]))
                else:
                    print(ind+var+u' = '+avar + u" =" + G2mth.ValEsd(val,0))
        # create a plot for each parametric variable
        for fitnum,obj in enumerate(eqObjList):
            calcobj = G2obj.ExpressionCalcObj(obj)
            # lookup dependent var position
            indx = colLabels.index(obj.GetDepVar())
            # assemble a list of the independent variables
            indepVars = obj.GetIndependentVars()
            # loop over each datapoint
            fitvals = []
            for j,row in enumerate(zip(*G2frame.colList)):
                calcobj.SetupCalc({var:row[i] for i,var in enumerate(colLabels) if var in indepVars})
                fitvals.append(calcobj.EvalExpression())
            G2plt.PlotSelectedSequence(G2frame,[indx],GetColumnInfo,SelectXaxis,fitnum,fitvals)

    def SingleParEqFit(eqObj):
        DoParEqFit(None,eqObj)

    def DelParFitEq(event):
        'Ask the user to select function to delete'
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in data['SeqParFitEqList']]
        selected = G2G.ItemSelector(
            txtlst,G2frame,
            multiple=True,
            title='Select a parametric equation(s) to remove',
            header='Delete equation')
        if selected is None: return
        data['SeqParFitEqList'] = [obj for i,obj in enumerate(data['SeqParFitEqList']) if i not in selected]
        EnableParFitEqMenus()
        if data['SeqParFitEqList']: DoParEqFit(event)

    def EditParFitEq(event):
        'Edit an existing parametric equation'
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in data['SeqParFitEqList']]
        if len(txtlst) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                txtlst,G2frame,
                multiple=False,
                title='Select a parametric equation to edit',
                header='Edit equation')
        if selected is not None:
            dlg = G2exG.ExpressionDialog(G2frame.dataDisplay,VarDict,
                data['SeqParFitEqList'][selected],depVarDict=VarDict,
                header="Edit the formula for this minimization function",
                ExtraButton=['Fit',SingleParEqFit],wildCard=False)
            newobj = dlg.Show(True)
            if newobj:
                data['SeqParFitEqList'][selected] = newobj
                EnableParFitEqMenus()
            if data['SeqParFitEqList']: DoParEqFit(event)

    def AddNewParFitEq(event):
        'Create a new parametric equation to be fit to sequential results'

        # compile the variable names used in previous freevars to avoid accidental name collisions
        usedvarlist = []
        for obj in data['SeqParFitEqList']:
            for var in obj.freeVars:
                if obj.freeVars[var][0] not in usedvarlist: usedvarlist.append(obj.freeVars[var][0])

        dlg = G2exG.ExpressionDialog(G2frame.dataDisplay,VarDict,depVarDict=VarDict,
            header='Define an equation to minimize in the parametric fit',
            ExtraButton=['Fit',SingleParEqFit],usedVars=usedvarlist,wildCard=False)
        obj = dlg.Show(True)
        dlg.Destroy()
        if obj:
            data['SeqParFitEqList'].append(obj)
            EnableParFitEqMenus()
            if data['SeqParFitEqList']: DoParEqFit(event)

    def CopyParFitEq(event):
        'Copy an existing parametric equation to be fit to sequential results'
        # compile the variable names used in previous freevars to avoid accidental name collisions
        usedvarlist = []
        for obj in data['SeqParFitEqList']:
            for var in obj.freeVars:
                if obj.freeVars[var][0] not in usedvarlist: usedvarlist.append(obj.freeVars[var][0])
        txtlst = [obj.GetDepVar()+' = '+obj.expression for obj in data['SeqParFitEqList']]
        if len(txtlst) == 1:
            selected = 0
        else:
            selected = G2G.ItemSelector(
                txtlst,G2frame,
                multiple=False,
                title='Select a parametric equation to copy',
                header='Copy equation')
        if selected is not None:
            newEqn = copy.deepcopy(data['SeqParFitEqList'][selected])
            for var in newEqn.freeVars:
                newEqn.freeVars[var][0] = G2obj.MakeUniqueLabel(newEqn.freeVars[var][0],usedvarlist)
            dlg = G2exG.ExpressionDialog(
                G2frame.dataDisplay,VarDict,newEqn,depVarDict=VarDict,
                header="Edit the formula for this minimization function",
                ExtraButton=['Fit',SingleParEqFit],wildCard=False)
            newobj = dlg.Show(True)
            if newobj:
                data['SeqParFitEqList'].append(newobj)
                EnableParFitEqMenus()
            if data['SeqParFitEqList']: DoParEqFit(event)

    def GridSetToolTip(row,col):
        '''Routine to show standard uncertainties for each element in table
        as a tooltip
        '''
        if G2frame.colSigs[col]:
            if G2frame.colSigs[col][row] == -0.1: return 'frozen'
            return u'\u03c3 = '+str(G2frame.colSigs[col][row])
        return ''

    def GridColLblToolTip(col):
        '''Define a tooltip for a column. This will be the user-entered value
        (from data['variableLabels']) or the default name
        '''
        if col < 0 or col > len(colLabels):
            print ('Illegal column #%d'%col)
            return
        var = colLabels[col]
        return variableLabels.get(var,G2obj.fmtVarDescr(var))

    def DoSequentialExport(event):
        '''Event handler for all Sequential Export menu items
        '''
        if event.GetId() == G2G.wxID_XPORTSEQFCIF:
            G2IO.ExportSequentialFullCIF(G2frame,data,Controls)
            return
        vals = G2frame.dataWindow.SeqExportLookup.get(event.GetId())
        if vals is None:
            print('Error: Id not found. This should not happen!')
            return
        G2IO.ExportSequential(G2frame,data,*vals)

    def onSelectSeqVars(event):
        '''Select which variables will be shown in table'''
        hides = [saveColLabels[2:].index(item) for item in G2frame.SeqTblHideList if
                     item in saveColLabels[2:]]
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select columns to hide',
                'Hide columns',saveColLabels[2:])
        dlg.SetSelections(hides)
        if dlg.ShowModal() == wx.ID_OK:
            G2frame.SeqTblHideList = [saveColLabels[2:][sel] for sel in dlg.GetSelections()]
            dlg.Destroy()
            UpdateSeqResults(G2frame,data,G2frame.dataDisplay.GetSize()) # redisplay variables
        else:
            dlg.Destroy()

    def OnCellChange(event):
        c = event.GetCol()
        if c != 1: return
        r = event.GetRow()
        val = G2frame.SeqTable.GetValue(r,c)
        data['Use'][r] = val
        G2frame.SeqTable.SetValue(r,c, val)

    def OnSelectUpdate(event):
        '''Update all phase parameters from a selected column in the Sequential Table.
        If no histogram is selected (or more than one), ask the user to make a selection.

        Loosely based on :func:`GSASIIstrIO.SetPhaseData`
        #TODO effect of ISO modes?
        '''
        rows = G2frame.dataDisplay.GetSelectedRows()
        if len(rows) == 1:
            sel = rows[0]
        else:
            dlg = G2G.G2SingleChoiceDialog(G2frame, 'Select a histogram to\nupdate phase from',
                                           'Select row',histNames)
            if dlg.ShowModal() == wx.ID_OK:
                sel = dlg.GetSelection()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return
        parmDict = data[histNames[sel]]['parmDict']
        Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
        for phase in Phases:
            print('Updating {} from Seq. Ref. row {}'.format(phase,histNames[sel]))
            Phase = Phases[phase]
            General = Phase['General']
            SGData = General['SGData']
            Atoms = Phase['Atoms']
            cell = General['Cell']
            pId = Phase['pId']
            pfx = str(pId)+'::'
            # there should not be any changes to the cell because those terms are not refined
            A,sigA = G2stIO.cellFill(pfx,SGData,parmDict,{})
            cell[1:7] = G2lat.A2cell(A)
            cell[7] = G2lat.calc_V(A)
            textureData = General['SH Texture']
            if textureData['Order']:
                for name in ['omega','chi','phi']:
                    aname = pfx+'SH '+name
                    textureData['Sample '+name][1] = parmDict[aname]
                for name in textureData['SH Coeff'][1]:
                    aname = pfx+name
                    textureData['SH Coeff'][1][name] = parmDict[aname]
            ik = 6  #for Pawley stuff below
            if General.get('Modulated',False):
                ik = 7
            # how are these updated?
            #General['SuperVec']
            #RBModels = Phase['RBModels']
            if Phase['General'].get('doPawley'):
                pawleyRef = Phase['Pawley ref']
                for i,refl in enumerate(pawleyRef):
                    key = pfx+'PWLref:'+str(i)
                    refl[ik] = parmDict[key]
#                    if key in sigDict:  #TODO: error here sigDict not defined. What was intended?
#                        refl[ik+1] = sigDict[key]
#                    else:
#                        refl[ik+1] = 0
                continue
            cx,ct,cs,cia = General['AtomPtrs']
            for i,at in enumerate(Atoms):
                names = {cx:pfx+'Ax:'+str(i),cx+1:pfx+'Ay:'+str(i),cx+2:pfx+'Az:'+str(i),cx+3:pfx+'Afrac:'+str(i),
                    cia+1:pfx+'AUiso:'+str(i),cia+2:pfx+'AU11:'+str(i),cia+3:pfx+'AU22:'+str(i),cia+4:pfx+'AU33:'+str(i),
                    cia+5:pfx+'AU12:'+str(i),cia+6:pfx+'AU13:'+str(i),cia+7:pfx+'AU23:'+str(i),
                    cx+4:pfx+'AMx:'+str(i),cx+5:pfx+'AMy:'+str(i),cx+6:pfx+'AMz:'+str(i)}
                for ind in range(cx,cx+4):
                    at[ind] = parmDict[names[ind]]
                if at[cia] == 'I':
                    at[cia+1] = parmDict[names[cia+1]]
                else:
                    for ind in range(cia+2,cia+8):
                        at[ind] = parmDict[names[ind]]
                if General['Type'] == 'magnetic':
                    for ind in range(cx+4,cx+7):
                        at[ind] = parmDict[names[ind]]
                ind = General['AtomTypes'].index(at[ct])
                if General.get('Modulated',False):
                    AtomSS = at[-1]['SS1']
                    waveType = AtomSS['waveType']
                    for Stype in ['Sfrac','Spos','Sadp','Smag']:
                        Waves = AtomSS[Stype]
                        for iw,wave in enumerate(Waves):
                            stiw = str(i)+':'+str(iw)
                            if Stype == 'Spos':
                                if waveType in ['ZigZag','Block',] and not iw:
                                    names = ['Tmin:'+stiw,'Tmax:'+stiw,'Xmax:'+stiw,'Ymax:'+stiw,'Zmax:'+stiw]
                                else:
                                    names = ['Xsin:'+stiw,'Ysin:'+stiw,'Zsin:'+stiw,
                                        'Xcos:'+stiw,'Ycos:'+stiw,'Zcos:'+stiw]
                            elif Stype == 'Sadp':
                                names = ['U11sin:'+stiw,'U22sin:'+stiw,'U33sin:'+stiw,
                                    'U12sin:'+stiw,'U13sin:'+stiw,'U23sin:'+stiw,
                                    'U11cos:'+stiw,'U22cos:'+stiw,'U33cos:'+stiw,
                                    'U12cos:'+stiw,'U13cos:'+stiw,'U23cos:'+stiw]
                            elif Stype == 'Sfrac':
                                if 'Crenel' in waveType and not iw:
                                    names = ['Fzero:'+stiw,'Fwid:'+stiw]
                                else:
                                    names = ['Fsin:'+stiw,'Fcos:'+stiw]
                            elif Stype == 'Smag':
                                names = ['MXsin:'+stiw,'MYsin:'+stiw,'MZsin:'+stiw,
                                    'MXcos:'+stiw,'MYcos:'+stiw,'MZcos:'+stiw]
                            for iname,name in enumerate(names):
                                AtomSS[Stype][iw][0][iname] = parmDict[pfx+name]

    def OnEditSelectPhaseVars(event):
        '''Select phase parameters in a selected histogram in a sequential
        fit. This allows the user to set their value(s)
        '''
        rows = G2frame.dataDisplay.GetSelectedRows()
        if len(rows) >= 1:
            selRows = rows
        else:
            dlg = G2G.G2MultiChoiceDialog(G2frame, 'Select histogram(s) to update\nphase parameters',
                                           'Select row',histNames)
            if dlg.ShowModal() == wx.ID_OK:
                selRows = dlg.GetSelections()
            else:
                selRows = []
            dlg.Destroy()
            if len(selRows) == 0: return
        parmDict = data[histNames[selRows[0]]]['parmDict']
        # narrow down to items w/o a histogram & having float values
        phaseKeys = [i for i in parmDict if ':' in i and i.split(':')[1] == '']
        phaseKeys = [i for i in phaseKeys if type(parmDict[i]) not in (int,str,bool)]
        if len(selRows) == 1:
            lbl = "\nin {}      ".format(histNames[selRows[0]])
        else:
            lbl = "\nin {} histograms".format(len(selRows))
        dlg = G2G.G2MultiChoiceDialog(G2frame, 'Choose phase parmDict item(s) to set'+lbl,
                                      'Choose items to edit', phaseKeys)
        if dlg.ShowModal() == wx.ID_OK:
            select = dlg.GetSelections()
            dlg.Destroy()
        else:
            dlg.Destroy()
            return
        if len(select) == 0: return
        l = [phaseKeys[i] for i in select]
        d = {i:parmDict[i] for i in l}
        val = G2G.CallScrolledMultiEditor(G2frame,len(l)*[d],l,l,CopyButton=True)
        if val:
            for sel in selRows:
                parmDict = data[histNames[sel]]['parmDict']
                for key in d: # update values shown in table
                    if parmDict[key] == d[key]: continue
                    if key in data[histNames[sel]]['varyList']:
                        i = data[histNames[sel]]['varyList'].index(key)
                        data[histNames[sel]]['variables'][i] = d[key]
                        data[histNames[sel]]['sig'][i] = 0
                    if key in data[histNames[sel]].get('depParmDict',{}):
                        data[histNames[sel]]['depParmDict'][key] = (d[key],0)
                parmDict.update(d) # update values used in next fit
        wx.CallAfter(UpdateSeqResults,G2frame,data) # redisplay variables
        return

#---- UpdateSeqResults: start processing sequential results here ##########
    # lookup table for unique cell parameters by symmetry
    if not data:
        print ('No sequential refinement results')
        return
    variableLabels = data.get('variableLabels',{})
    data['variableLabels'] = variableLabels
    Histograms,Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
    ifPWDR = True
    if not len(Histograms) and not len(Phases):   #SASD, REFD, PDF or IMG histogrms not PWDR or HKLF
        histNames = [name for name in data['histNames']]
        Controls = {}
        ifPWDR = False
    else:
        Controls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame,G2frame.root,'Controls'))
        # create a place to store Pseudo Vars & Parametric Fit functions, if not present
        if 'SeqPseudoVars' not in data: data['SeqPseudoVars'] = {}
        if 'SeqParFitEqList' not in data: data['SeqParFitEqList'] = []
        histNames = [name for name in data['histNames'] if name in data]
    if G2frame.dataDisplay:
        G2frame.dataDisplay.Destroy()
    G2frame.GetStatusBar().SetStatusText("Select column to export; LMB/RMB column to plot data/change label; LMB/RMB on row for PWDR/Covariance plot",1)
    sampleParms = GetSampleParms()

    # get unit cell & symmetry for all phases & initial stuff for later use
    RecpCellTerms = {}
    SGdata = {}
    uniqCellIndx = {}
    initialCell = {}
    RcellLbls = {}
    zeroDict = {}
    for phase in Phases:
        phasedict = Phases[phase]
        pId = phasedict['pId']
        pfx = str(pId)+'::' # prefix for A values from phase
        RcellLbls[pId] = [pfx+'A'+str(i) for i in range(6)]
        RecpCellTerms[pId] = G2lat.cell2A(phasedict['General']['Cell'][1:7])
        zeroDict[pId] = dict(zip(RcellLbls[pId],6*[0.,]))
        SGdata[pId] = phasedict['General']['SGData']
        laue = SGdata[pId]['SGLaue']
        if laue == '2/m':
            laue += SGdata[pId]['SGUniq']
        for symlist,celllist in G2lat.UniqueCellByLaue:
            if laue in symlist:
                uniqCellIndx[pId] = celllist
                break
        else: # should not happen
            uniqCellIndx[pId] = list(range(6))
        for i in uniqCellIndx[pId]:
            initialCell[str(pId)+'::A'+str(i)] =  RecpCellTerms[pId][i]

    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.PWDRMenu)
    G2gd.SetDataMenuBar(G2frame,G2frame.dataWindow.SequentialMenu)
    G2frame.Bind(wx.EVT_MENU, OnSelectUse, id=G2G.wxID_SELECTUSE)
    G2frame.Bind(wx.EVT_MENU, OnRenameSelSeq, id=G2G.wxID_RENAMESEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnSaveSelSeq, id=G2G.wxID_SAVESEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnSaveSelSeqCSV, id=G2G.wxID_SAVESEQSELCSV)
    G2frame.Bind(wx.EVT_MENU, OnSaveSeqCSV, id=G2G.wxID_SAVESEQCSV)
    G2frame.Bind(wx.EVT_MENU, OnPlotSelSeq, id=G2G.wxID_PLOTSEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnAveSelSeq, id=G2G.wxID_AVESEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnSelectUpdate, id=G2G.wxID_UPDATESEQSEL)
    G2frame.Bind(wx.EVT_MENU, OnEditSelectPhaseVars, id=G2G.wxID_EDITSEQSELPHASE)
    G2frame.Bind(wx.EVT_MENU, onSelectSeqVars, id=G2G.wxID_ORGSEQINC)
    G2frame.Bind(wx.EVT_MENU, AddNewPseudoVar, id=G2G.wxID_ADDSEQVAR)
    G2frame.Bind(wx.EVT_MENU, AddNewDistPseudoVar, id=G2G.wxID_ADDSEQDIST)
    G2frame.Bind(wx.EVT_MENU, AddNewAnglePseudoVar, id=G2G.wxID_ADDSEQANGLE)
    G2frame.Bind(wx.EVT_MENU, DelPseudoVar, id=G2G.wxID_DELSEQVAR)
    G2frame.Bind(wx.EVT_MENU, EditPseudoVar, id=G2G.wxID_EDITSEQVAR)
    G2frame.Bind(wx.EVT_MENU, AddNewParFitEq, id=G2G.wxID_ADDPARFIT)
    G2frame.Bind(wx.EVT_MENU, CopyParFitEq, id=G2G.wxID_COPYPARFIT)
    G2frame.Bind(wx.EVT_MENU, DelParFitEq, id=G2G.wxID_DELPARFIT)
    G2frame.Bind(wx.EVT_MENU, EditParFitEq, id=G2G.wxID_EDITPARFIT)
    G2frame.Bind(wx.EVT_MENU, DoParEqFit, id=G2G.wxID_DOPARFIT)

    for id in G2frame.dataWindow.SeqExportLookup:
        G2frame.Bind(wx.EVT_MENU, DoSequentialExport, id=id)
    G2frame.Bind(wx.EVT_MENU, OnSaveSeqCSV, id=G2G.wxID_XPORTSEQCSV)
    G2frame.Bind(wx.EVT_MENU, DoSequentialExport, id=G2G.wxID_XPORTSEQFCIF)

    EnablePseudoVarMenus()
    EnableParFitEqMenus()

    combinedVaryList = []  # list of all varied parameters ; used for column headings
    atomsVaryList = {}     # dict of atom coords varied in any histogram, includes dependent params
                           # key is atom param name, value is esd parm name
    firstValueDict = {}    # first value for each parameter; used to create VarDict for parametric fit pseudo vars GUI
    foundHistNames = []    # histograms to be used in sequential table
    maxPWL = 5             # number of Pawley vars to show
    missing = 0
    newCellDict = {}
    PSvarDict = {} # contains 1st value for each parameter in parmDict
                   # needed for creating & editing pseudovars
    PSvarDict.update(sampleParms)

    # scan through histograms to see what are available and to make a
    # list of all varied parameters; also create a dict that has the
    for i,name in enumerate(histNames):
        if name not in data:
            if missing < 5:
                print(" Warning: "+name+" not found")
            elif missing == 5:
                print (' Warning: more are missing')
            missing += 1
            continue
        foundHistNames.append(name)
        for var,val,sig in zip(data[name]['varyList'],data[name]['variables'],data[name]['sig']):
            svar = striphist(var,'*') # wild-carded
            if 'PWL' in svar:
                if int(svar.split(':')[-1]) > maxPWL:
                    continue
            if svar not in combinedVaryList:
                # add variables to list as they appear
                combinedVaryList.append(svar)
                firstValueDict[svar] = (val,sig)
            if '::dA' in var:
                atomsVaryList[var.replace('::dA','::A')] = var
                atomsVaryList.update({i.replace('::dA','::A'):i for i in data[name]['depParmDict'] if '::dA' in i})
        # get all refined cell terms
        newCellDict.update(data[name].get('newCellDict',{})) # N.B. These Dij vars are missing a histogram #
        # make sure 1st reference to each parm is in PseudoVar dict
        tmp = copy.deepcopy(data[name].get('parmDict',{}))
        tmp = {striphist(var,'*'):tmp[var] for var in tmp}  # replace histogram #s with "*"
        tmp.update(PSvarDict)
        PSvarDict = tmp

    if missing:
        print (' Warning: Total of %d data sets missing from sequential results'%(missing))

    # make dict of varied cell parameters equivalents
    ESDlookup = {} # provides the Dij term for each Ak term (where terms are refined)
    Dlookup = {} # provides the Ak term for each Dij term (where terms are refined)
    cellAlist = []
    for item in newCellDict:
        cellAlist.append(newCellDict[item][0])
        if item in data.get('varyList',[]):
            ESDlookup[newCellDict[item][0]] = item
            Dlookup[item] = newCellDict[item][0]
    # add coordinate equivalents to lookup table
    for parm in atomsVaryList:
        Dlookup[atomsVaryList[parm]] = parm
        ESDlookup[parm] = atomsVaryList[parm]
    combinedVaryList.sort()
    histNames = foundHistNames
    # prevVaryList = []
    posdict = {}    # defines position for each entry in table; for inner
                    # dict key is column number & value is parameter name
    histNumList = []
    for i,name in enumerate(histNames):
        if name in Histograms:
            histNumList.append(list(Histograms.keys()).index(name))
        # if prevVaryList != data[name]['varyList']: # this refinement has a different refinement list from previous
        #     prevVaryList = data[name]['varyList']
        posdict[name] = {}
        for var in data[name]['varyList']:
            svar = striphist(var,'*')
            if 'PWL' in svar:
                if int(svar.split(':')[-1]) > maxPWL:
                    continue
            posdict[name][combinedVaryList.index(svar)] = svar
    ####--  build up the data table by columns -----------------------------------------------
    nRows = len(histNames)
    if len(histNumList) != nRows:
        G2frame.colList = [list(range(nRows))]
    else:
        G2frame.colList = [histNumList]
    if len(data.get('Use',[])) != nRows:
        data['Use'] = nRows*[True]
    G2frame.colList += [data['Use']]
    G2frame.colSigs = [None,None,]
    colLabels = ['No.','Use',]
    Types = [wg.GRID_VALUE_LONG,wg.GRID_VALUE_BOOL,]
    # start with Rwp values
    if 'IMG ' not in histNames[0][:4]:
        G2frame.colList += [[data[name]['Rvals']['Rwp'] for name in histNames]]
        G2frame.colSigs += [None]
        colLabels += ['Rwp']
        Types += [wg.GRID_VALUE_FLOAT+':10,3',]
    if histNames[0][:4] not in ['SASD','IMG ','REFD','PDF ']:
        G2frame.colList += [[data[name]['Rvals']['GOF'] for name in histNames]]
        G2frame.colSigs += [None]
        colLabels += ['GOF']
        Types += [wg.GRID_VALUE_FLOAT+':10,3',]
    # add % change in Chi^2 in last cycle
    if histNames[0][:4] not in ['SASD','IMG ','REFD'] and Controls.get('ShowCell'):
        G2frame.colList += [[100.*data[name]['Rvals'].get('DelChi2',-1) for name in histNames]]
        G2frame.colSigs += [None]
        colLabels += [u'\u0394\u03C7\u00B2 (%)']
        Types += [wg.GRID_VALUE_FLOAT+':10,5',]
    deltaChiCol = len(colLabels)-1
    # frozen variables?
    if 'parmFrozen' in Controls:
        f = [len(Controls['parmFrozen'].get(h,[])) for h in histNames]
        if any(f):
            G2frame.colList += [f]
            G2frame.colSigs += [None]
            colLabels += ['frozen']
            Types += [wg.GRID_VALUE_LONG]
    # add changing sample parameters to table
    for key in sampleParms:
        G2frame.colList += [sampleParms[key]]
        G2frame.colSigs += [None]
        colLabels += [key]
        Types += [wg.GRID_VALUE_FLOAT,]
    sampleDict = {}
    for i,name in enumerate(histNames):
        sampleDict[name] = dict(zip(sampleParms.keys(),[sampleParms[key][i] for key in sampleParms.keys()]))
    # add unique cell parameters
    if Controls.get('ShowCell',False) and len(newCellDict):
        phaseLookup = {Phases[phase]['pId']:phase for phase in Phases}
        for pId in sorted(RecpCellTerms):
            pfx = str(pId)+'::' # prefix for A values from phase
            cells = []
            cellESDs = []
            colLabels += [pfx+G2lat.cellUlbl[i] for i in uniqCellIndx[pId]]
            colLabels += [pfx+'Vol']
            Types += (len(uniqCellIndx[pId]))*[wg.GRID_VALUE_FLOAT+':10,5',]
            Types += [wg.GRID_VALUE_FLOAT+':10,3',]
            Albls = [pfx+'A'+str(i) for i in range(6)]
            for name in histNames:
                if name not in Histograms: continue
                hId = Histograms[name]['hId']
                phfx = '%d:%d:'%(pId,hId)
                esdLookUp = {}
                dLookup = {}
                for item in data[name]['newCellDict']:
                    if phfx+item.split('::')[1] in data[name]['varyList']:
                        esdLookUp[newCellDict[item][0]] = item
                        dLookup[item] = newCellDict[item][0]
                covData = {'varyList': [dLookup.get(striphist(v),v) for v in data[name]['varyList']],
                    'covMatrix': data[name]['covMatrix']}
                A = RecpCellTerms[pId][:] # make copy of starting A values
                # update with refined values
                for i,j in enumerate(('D11','D22','D33','D12','D13','D23')):
                    var = str(pId)+'::A'+str(i)
                    Dvar = str(pId)+':'+str(hId)+':'+j
                    # apply Dij value if non-zero
                    if Dvar in data[name]['parmDict']:
                        if data[name]['parmDict'][Dvar] != 0.0:
                            A[i] += data[name]['parmDict'][Dvar]
                    # override with fit result if is Dij varied
                    if var in cellAlist:
                        try:
                            A[i] = data[name]['newCellDict'][esdLookUp[var]][1] # get refined value
                        except KeyError:
                            pass
                # apply symmetry
                cellDict = dict(zip(Albls,A))
                try:    # convert to direct cell
                    A,zeros = G2stIO.cellFill(pfx,SGdata[pId],cellDict,zeroDict[pId])
                    c = G2lat.A2cell(A)
                    vol = G2lat.calc_V(A)
                    cE = G2lat.getCellEsd(pfx,SGdata[pId],A,covData)
                except:
                    c = 6*[None]
                    cE = 6*[None]
                    vol = None
                # add only unique values to table
                if name in Phases[phaseLookup[pId]]['Histograms']:
                    cells += [[c[i] for i in uniqCellIndx[pId]]+[vol]]
                    cellESDs += [[cE[i] for i in uniqCellIndx[pId]]+[cE[-1]]]
                    # add in direct cell terms to PseudoVar dict
                    tmp = dict(zip([pfx+G2lat.cellUlbl[i] for i in uniqCellIndx[pId]]+[pfx+'Vol'],
                                     [c[i] for i in uniqCellIndx[pId]]+[vol]))
                    tmp.update(PSvarDict)
                    PSvarDict = tmp
                else:
                    cells += [[None for i in uniqCellIndx[pId]]+[None]]
                    cellESDs += [[None for i in uniqCellIndx[pId]]+[None]]
            G2frame.colList += zip(*cells)
            G2frame.colSigs += zip(*cellESDs)

    # get ISODISTORT labels
    ISOlist = []
    for phase in Phases:
        ISOlist += [i.varname() for i in Phases[phase].get('ISODISTORT',{}).get('G2ModeList',[])
                       if i.varname() not in ISOlist]
    # set labels for columns of data table
    ISOcols = {}  # ISODISTORT modes
    for i,lbl in enumerate(combinedVaryList):
        if 'nv-' in lbl:
            nlbl = lbl.replace('::nv-','::')
            if nlbl in ISOlist:
                lbl = nlbl
                ISOcols[lbl] = i
        colLabels.append(lbl)
    Types += len(combinedVaryList)*[wg.GRID_VALUE_FLOAT,]
    vals = []
    esds = []
    varsellist = None        # will be a list of variable names in the order they are selected to appear
    # tabulate values for each hist, leaving None for blank columns
    for ih,name in enumerate(histNames):
        varsellist = [posdict[name].get(i) for i in range(len(combinedVaryList))]
        # translate variable names to how they will be used in the headings
        vs = [striphist(v,'*') for v in data[name]['varyList']]
        # determine the index for each column (or None) in the data[]['variables'] and ['sig'] lists
        sellist = [vs.index(v) if v is not None else None for v in varsellist]
        #sellist = [i if striphist(v,'*') in varsellist else None for i,v in enumerate(data[name]['varyList'])]
        if not varsellist: raise Exception()
        vals.append([data[name]['variables'][s] if s is not None else None for s in sellist])
        #replace mode displacement shift with value; esd applies to both
        for pname in ISOcols:
            if pname in data[name]['parmDict']:
                vals[ih][ISOcols[pname]] = data[name]['parmDict'][pname]
        esds.append([data[name]['sig'][s] if s is not None else None for s in sellist])
    G2frame.colList += zip(*vals)
    G2frame.colSigs += zip(*esds)

    # add refined atom parameters to table
    for parm in sorted(atomsVaryList):
        vals = []
        sigs = []
        aprm = atomsVaryList[parm]
        for name in histNames:
            if aprm in data[name]['varyList']:
                i = data[name]['parmDict'][parm]
                j = data[name]['sig'][data[name]['varyList'].index(aprm)]
            elif aprm in data[name]['depParmDict']:
                i = data[name]['parmDict'][parm]
                j = data[name]['depParmDict'][aprm][1]
            else:
                i = j = None
            vals.append(i)
            sigs.append(j)
        colLabels.append(parm)
        Types += [wg.GRID_VALUE_FLOAT+':10,5',]
        G2frame.colSigs += [sigs]
        G2frame.colList += [vals]

    # tabulate dependent parameters from constraints, removing histogram numbers from
    # parm label, if needed. Add the dependent variables to the table
    depValDict = {}
    for name in histNames:
        for var in data[name].get('depParmDict',{}):
            if '::dA' in var: continue
            val,sig = data[name]['depParmDict'][var]
            svar = striphist(var,'*')
            if svar not in depValDict:
               depValDict[svar] = {}
            depValDict[svar][name] = (val,sig)
    for svar in sorted(depValDict):
        vals = []
        sigs = []
        for name in histNames:
            if name in depValDict[svar]:
                i,j = depValDict[svar][name]
            else:
                i = j = None
            vals.append(i)
            sigs.append(j)
        colLabels.append(svar)
        Types += [wg.GRID_VALUE_FLOAT+':10,5',]
        G2frame.colSigs += [sigs]
        G2frame.colList += [vals]

    # evaluate Pseudovars, their ESDs and add them to grid
    # this should be reworked so that the eval dict is created once and as noted below
    for expr in data['SeqPseudoVars']:
        obj = data['SeqPseudoVars'][expr]
        calcobj = G2obj.ExpressionCalcObj(obj)
        valList = []
        esdList = []
        for seqnum,name in enumerate(histNames):
            sigs = data[name]['sig']
            G2mv.InitVars()
#            parmDict = data[name].get('parmDict')
            parmDict = data[name]['parmDict']
            constraintInfo = data[name].get('constraintInfo',[[],[],{},[],seqnum])
            groups,parmlist,constrDict,fixedList,ihst = constraintInfo
            varyList = data[name]['varyList']
            msg = G2mv.EvaluateMultipliers(constrDict,parmDict)
            if msg:
                print('Unable to interpret multiplier(s) for',name,':',msg)
                continue
            G2mv.GenerateConstraints(varyList,constrDict,fixedList,parmDict,
                                     seqHistNum=ihst,raiseException=False)
            if 'Dist' in expr:
                derivs = G2mth.CalcDistDeriv(obj.distance_dict,obj.distance_atoms, parmDict)
                pId = obj.distance_dict['pId']
                aId,bId = obj.distance_atoms
                varyNames = ['%d::dA%s:%d'%(pId,ip,aId) for ip in ['x','y','z']]
                varyNames += ['%d::dA%s:%d'%(pId,ip,bId) for ip in ['x','y','z']]
                VCoV = G2mth.getVCov(varyNames,varyList,data[name]['covMatrix'])
                esdList.append(np.sqrt(np.inner(derivs,np.inner(VCoV,derivs.T)) ))
            elif 'Angle' in expr:
                derivs = G2mth.CalcAngleDeriv(obj.angle_dict,obj.angle_atoms, parmDict)
                pId = obj.angle_dict['pId']
                aId,bId = obj.angle_atoms
                varyNames = ['%d::dA%s:%d'%(pId,ip,aId) for ip in ['x','y','z']]
                varyNames += ['%d::dA%s:%d'%(pId,ip,bId[0]) for ip in ['x','y','z']]
                varyNames += ['%d::dA%s:%d'%(pId,ip,bId[1]) for ip in ['x','y','z']]
                VCoV = G2mth.getVCov(varyNames,varyList,data[name]['covMatrix'])
                esdList.append(np.sqrt(np.inner(derivs,np.inner(VCoV,derivs.T)) ))
            else:
                derivs = np.array(  # TODO: this needs to be reworked
                    [EvalPSvarDeriv(calcobj,parmDict.copy(),sampleDict[name],var,ESD)
                     for var,ESD in zip(varyList,sigs)])
                # needs work: calls calcobj.SetupCalc each call time
                # integrate into G2obj.ExpressionCalcObj
                if None in list(derivs):
                    esdList.append(None)
                else:
                    esdList.append(np.sqrt(
                        np.inner(derivs,np.inner(data[name]['covMatrix'],derivs.T)) ))
            psDict = parmDict.copy()
            psDict.update(sampleDict[name])
            try:
                UpdateParmDict(psDict)
            except:
                print('UpdateParmDict error on histogram',name)
            calcobj.UpdateDict(psDict)
            valList.append(calcobj.EvalExpression())
#            if calcobj.su is not None: esdList[-1] = calcobj.su
        if not esdList:
            esdList = None
        G2frame.colList += [valList]
        G2frame.colSigs += [esdList]
        colLabels += [expr]
        Types += [wg.GRID_VALUE_FLOAT+':10,5']
    #---- table build done -------------------------------------------------------------

    # clean up the PseudoVars dict by reomving dA[xyz] & Dij
    remDij =   re.compile('[0-9]+:[0-9]*:D[123][123]')
    remdAxyz = re.compile('[0-9]+::dA[xyz]:[0-9]+')
    PSvarDict = {i:PSvarDict[i] for i in PSvarDict if not (remDij.match(i) or remdAxyz.match(i))}

    # remove selected columns from table
    saveColLabels = colLabels[:]
    if G2frame.SeqTblHideList is None:      #set default hides
        G2frame.SeqTblHideList = [item for item in saveColLabels if 'Back' in item and ifPWDR]
        G2frame.SeqTblHideList += [item for item in saveColLabels if 'dA' in item]
        G2frame.SeqTblHideList += [item for item in saveColLabels if ':*:D' in item]
    #******************************************************************************
    # create a set of values for example evaluation of parametric equation fitting
    VarDict = {}
    for i,var in enumerate(colLabels):
        if var in ['Use','Rwp',u'\u0394\u03C7\u00B2 (%)']: continue
        if G2frame.colList[i][0] is None:
            val,sig = firstValueDict.get(var,[None,None])
        elif G2frame.colSigs[i]:
            val,sig = G2frame.colList[i][0],G2frame.colSigs[i][0]
        else:
            val,sig = G2frame.colList[i][0],None
        if striphist(var) not in Dlookup:
            VarDict[var] = val
    # add recip cell coeff. values
    VarDict.update({var:val for var,val in newCellDict.values()})

    # remove items to be hidden from table
    for l in reversed(range(len(colLabels))):
        if colLabels[l] in G2frame.SeqTblHideList:
            del colLabels[l]
            del Types[l]
            del G2frame.colList[l]
            del G2frame.colSigs[l]
            if deltaChiCol == l:
                deltaChiCol = None

    # make a copy of the column labels substituting alternate labels when defined
    displayLabels = colLabels[:]
    for i,l in enumerate(colLabels):
        if l in variableLabels:
            displayLabels[i] = variableLabels[l]

    G2frame.dataWindow.ClearData()
    G2frame.dataWindow.currentGrids = []
    G2frame.dataDisplay = G2G.GSGrid(parent=G2frame.dataWindow)
    G2frame.dataDisplay.SetScrollRate(10,10)
    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label='Sequential results:'),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add(G2frame.dataDisplay,1,wx.EXPAND,1)
    G2frame.dataWindow.SetSizer(mainSizer)
    if histNames[0].startswith('PWDR'):
        #rowLabels = [str(i)+': '+l[5:30] for i,l in enumerate(histNames)]
        rowLabels = [l[5:] for i,l in enumerate(histNames)]
    else:
        rowLabels = histNames
    G2frame.SeqTable = G2G.Table([list(cl) for cl in zip(*G2frame.colList)],     # convert from columns to rows
        colLabels=displayLabels,rowLabels=rowLabels,types=Types)
    G2frame.dataDisplay.SetTable(G2frame.SeqTable, True)
    # make all Use editable all others ReadOnly
    for c in range(len(colLabels)):
        for r in range(nRows):
            if c == 1:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=False)
            else:
                G2frame.dataDisplay.SetReadOnly(r,c,isReadOnly=True)
    if 'phoenix' in wx.version():
        G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGED, OnCellChange)
    else:
        G2frame.dataDisplay.Bind(wg.EVT_GRID_CELL_CHANGE, OnCellChange)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_LEFT_CLICK, PlotLeftSelect)
    G2frame.dataDisplay.Bind(wg.EVT_GRID_LABEL_RIGHT_CLICK, PlotRightSelect)
    G2frame.dataDisplay.SetRowLabelSize(8*len(histNames[0]))       #pretty arbitrary 8
    G2frame.dataDisplay.SetMargins(0,0)
    G2frame.dataDisplay.AutoSizeColumns(False)
    # highlight unconverged shifts
    if histNames[0][:4] not in ['SASD','IMG ','REFD',] and deltaChiCol is not None:
        for row,name in enumerate(histNames):
            if name not in Controls.get('Seq Data',{}):
                G2frame.dataDisplay.SetCellTextColour(row,0,wx.Colour(255,0,0))
            deltaChi = G2frame.SeqTable.GetValue(row,deltaChiCol)
            try:
                if deltaChi > 10.:
                    G2frame.dataDisplay.SetCellStyle(row,deltaChiCol,color=wx.Colour(255,0,0))
                elif deltaChi > 1.0:
                    G2frame.dataDisplay.SetCellStyle(row,deltaChiCol,color=wx.Colour(255,255,0))
            except:
                pass
    G2frame.dataDisplay.InstallGridToolTip(GridSetToolTip,GridColLblToolTip)
    #G2frame.dataDisplay.SendSizeEvent() # resize needed on mac
    #G2frame.dataDisplay.Refresh() # shows colored text on mac
    G2frame.dataWindow.SetDataSize()

###############################################################################################################
#UpdateClusterAnalysis: results
###############################################################################################################

def UpdateClusterAnalysis(G2frame,ClusData,shoNum=-1):
    import scipy.spatial.distance as SSD
    import scipy.cluster.hierarchy as SCH
    import scipy.cluster.vq as SCV
    try:
        import sklearn.cluster as SKC
        import sklearn.ensemble as SKE
        import sklearn.neighbors as SKN
        import sklearn.svm as SKVM
        import sklearn.metrics as SKM
        ClusData['SKLearn'] = True
    except:
        ClusData['SKLearn'] = False

    def FileSizer():

        def OnSelectData(event):

            def GetCaLimits(names):
                ''' scan through data selected for cluster analysis to find highest lower & lowest upper limits
                param: data dict: Cluster analysis info
                '''
                limits = [0.,1.e6]
                for name in names:
                    item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                    if 'PWDR' in name:
                        x = G2frame.GPXtree.GetItemPyData(item)[1][0]
                    else:
                        PDFControls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame, item,'PDF Controls'))
                        x = PDFControls['G(R)'][1][0]
                    limits = [max(np.min(x),limits[0]),min(np.max(x),limits[1])]
                return limits

            choices = G2gd.GetGPXtreeDataNames(G2frame,['PWDR','PDF '])
            if len(choices) == 0:
                G2G.G2MessageBox(G2frame,'No PWDR or PDF histograms found for cluster analysis.','No Histograms')
                return
            sel = []
            try:
                if 'Cluster Data' in ClusData:
                    sel = [choices.index(item) for item in ClusData['Cluster Data']['Files']]
            except ValueError:  #data changed somehow - start fresh
                sel = []
            dlg = G2G.G2MultiChoiceDialog(G2frame,
                'Select datasets to include.\n PWDR or PDF',
                'Cluster analysis data selection',choices)
            dlg.SetSelections(sel)
            names = []
            Type = ''
            if dlg.ShowModal() == wx.ID_OK:
                for sel in dlg.GetSelections():
                    if not Type:
                        Type = choices[sel].split()[0]
                    if Type != choices[sel].split()[0]:
                        G2G.G2MessageBox(G2frame,'Histogram types not all the same; revise selection','Histogram type mismatch')
                        return
                    names.append(choices[sel])
                ClusData['Files'] = names
                limits = GetCaLimits(names)
                ClusData['Limits'] = [copy.copy(limits),limits]
                ClusData['DataMatrix'] = []
                ClusData['ConDistMat'] = []
                ClusData['CLuZ'] = None
                ClusData['codes'] = None
                ClusData['plots'] = 'All'

            dlg.Destroy()
            G2frame.SetTitleByGPX()

            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        fileSizer = wx.BoxSizer(wx.HORIZONTAL)
        Type = 'PWDR'
        if len(ClusData['Files']):
            if 'PDF' in ClusData['Files'][0]:
                Type = 'PDF'
            lbl = 'Cluster Analysis with %d %s datasets: '%(len(ClusData['Files']),Type)
            ClusData['Type'] = Type
        else:
            lbl = 'No data selected for Cluster Analysis'
        fileSizer.Add(wx.StaticText(G2frame.dataWindow,label=lbl),0,WACV)
        selSeqData = wx.Button(G2frame.dataWindow,label=' Select datasets')
        selSeqData.Bind(wx.EVT_BUTTON,OnSelectData)
        fileSizer.Add(selSeqData,0,WACV)
        return fileSizer

    def LimitSizer():

        def CheckLimits(invalid,value,tc):
            #TODO this needs a check on ultimate size of data array; loop over names & count points?

            if ClusData['Limits'][1][1] < ClusData['Limits'][1][0]:
                ClusData['Limits'][1] = [ClusData['Limits'][1][1],ClusData['Limits'][1][0]]
            ClusData['DataMatrix'] = []
            ClusData['ConDistMat'] = []
            ClusData['CLuZ'] = None

            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        limitSizer = wx.BoxSizer(wx.HORIZONTAL)
        limitSizer.Add(wx.StaticText(G2frame.dataWindow,label='Enter cluster analysis data limits: '),0,WACV)
        limitSizer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,ClusData['Limits'][1],0,nDig=(10,3),
            xmin=ClusData['Limits'][0][0],xmax=ClusData['Limits'][0][1],OnLeave=CheckLimits),0,WACV)
        limitSizer.Add(G2G.ValidatedTxtCtrl(G2frame.dataWindow,ClusData['Limits'][1],1,nDig=(10,3),
            xmin=ClusData['Limits'][0][0],xmax=ClusData['Limits'][0][1],OnLeave=CheckLimits),0,WACV)
        return limitSizer

    def GetYMatSize():
        nData = 0
        for name in ClusData['Files']:
            item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
            if 'PWDR' in name:
                x = G2frame.GPXtree.GetItemPyData(item)[1][0]
            else:
                PDFControls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame, item,'PDF Controls'))
                x = PDFControls['G(R)'][1][0]
            iBeg = np.searchsorted(x,ClusData['Limits'][1][0])
            iFin = np.searchsorted(x,ClusData['Limits'][1][1])+1
            nData += (iFin-iBeg)
        return nData

    def OnMakeArray(event):
        Limits = ClusData['Limits'][1]
        Start = True
        nFiles = len(ClusData['Files'])
        CAmatrix = []
        try:
            for iname,name in enumerate(ClusData['Files']):
                item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
                if 'PWDR' in name:
                    x = G2frame.GPXtree.GetItemPyData(item)[1][0]
                    y = G2frame.GPXtree.GetItemPyData(item)[1][1]
                else:
                    PDFControls = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame, item,'PDF Controls'))
                    x = PDFControls['G(R)'][1][0]
                    y = PDFControls['G(R)'][1][1]
                iBeg = np.searchsorted(x,Limits[0])
                iFin = np.searchsorted(x,Limits[1])
                if Start:
                    CAmatrix = np.empty((nFiles,iFin-iBeg+1))
                    CAmatrix[iname] = y[iBeg:iFin+1]
                    Start = False
                else:
                    CAmatrix[iname] = y[iBeg:iFin+1]
        except ValueError:
            G2G.G2MessageBox(G2frame,
                'Data for %s is mismatched in length to those already processed or has different step size'%name,
                'No Cluster Analysis possible')
            return
        ClusData['DataMatrix'] = CAmatrix
        ClusData['ConDistMat'] = []
        ClusData['CLuZ'] = None
        wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

    def MethodSizer():

        def OnClusterMethod(event):
            ClusData['Method'] = method.GetValue()
            ClusData['ConDistMat'] = []
            ClusData['CLuZ'] = None
            OnCompute(event)
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        def OnCompute(event):
            if 'minkowski' in ClusData['Method']:
                ClusData['ConDistMat'] = SSD.pdist(ClusData['DataMatrix'],ClusData['Method'],p=int(ClusData['MinkP']))
            else:
                ClusData['ConDistMat'] = SSD.pdist(ClusData['DataMatrix'],ClusData['Method'])
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        def OnExponent(event):
            ClusData['MinkP'] = minp.GetValue()
            ClusData['ConDistMat'] = []
            ClusData['CLuZ'] = None
            OnCompute(event)
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        choice = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine',  \
            'euclidean', 'jensenshannon', 'minkowski', 'seuclidean',  'sqeuclidean']
        methsizer = wx.BoxSizer(wx.HORIZONTAL)
        methsizer.Add(wx.StaticText(G2frame.dataWindow,label='Select cluster analysis distance method: '),0,WACV)
        method = wx.ComboBox(G2frame.dataWindow,choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        method.SetValue(ClusData['Method'])
        method.Bind(wx.EVT_COMBOBOX, OnClusterMethod)
        methsizer.Add(method,0,WACV)
        if 'minkowski' in ClusData['Method']:
            methsizer.Add(wx.StaticText(G2frame.dataWindow,label=' exponent: '),0,WACV)
            choicep = ['1','2','3','4','10']
            minp = wx.ComboBox(G2frame.dataWindow,choices=choicep,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            minp.SetValue(ClusData['MinkP'])
            minp.Bind(wx.EVT_COMBOBOX, OnExponent)
            methsizer.Add(minp,0,WACV)
        compute = wx.Button(G2frame.dataWindow,label='Compute distance matrix')
        compute.Bind(wx.EVT_BUTTON,OnCompute)
        methsizer.Add(compute,0,WACV)
        return methsizer

    def HierSizer():

        def OnLinkMethod(event):
            ClusData['LinkMethod'] = method.GetValue()
            OnCompute(event)

        def OnOptOrd(event):
            ClusData['Opt Order'] = not ClusData['Opt Order']
            OnCompute(event)

        def OnCompute(event):
            ClusData['CLuZ'] = SCH.linkage(ClusData['ConDistMat'],method=ClusData['LinkMethod'],optimal_ordering=ClusData['Opt Order'])
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        hierSizer = wx.BoxSizer(wx.HORIZONTAL)
        hierSizer.Add(wx.StaticText(G2frame.dataWindow,label='Hierarchical clustering: Select linkage method: '),0,WACV)
        choice = ['single','complete','average','weighted','centroid','median','ward',]
        method = wx.ComboBox(G2frame.dataWindow,choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        method.SetValue(ClusData['LinkMethod'])
        method.Bind(wx.EVT_COMBOBOX, OnLinkMethod)
        hierSizer.Add(method,0,WACV)
        optOrd = wx.CheckBox(G2frame.dataWindow,label=' Optimal order? ')
        optOrd.Bind(wx.EVT_CHECKBOX,OnOptOrd)
        optOrd.SetValue(ClusData['Opt Order'])
        hierSizer.Add(optOrd,0,WACV)
        compute = wx.Button(G2frame.dataWindow,label='Compute')
        compute.Bind(wx.EVT_BUTTON,OnCompute)
        hierSizer.Add(compute,0,WACV)
        return hierSizer

    def kmeanSizer():

        def OnClusNum(event):
            ClusData['NumClust'] = int(numclust.GetValue())
            OnCompute(event)

        def OnCompute(event):
            whitMat = SCV.whiten(ClusData['DataMatrix'])
            codebook,dist = SCV.kmeans2(whitMat,ClusData['NumClust'])   #use K-means++
            ClusData['codes'],ClusData['dists'] = SCV.vq(whitMat,codebook)
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        kmeanssizer = wx.BoxSizer(wx.HORIZONTAL)
        choice = [str(i) for i in range(2,16)]
        kmeanssizer.Add(wx.StaticText(G2frame.dataWindow,label='K-means clustering: select number of clusters (2-15): '),0,WACV)
        numclust = wx.ComboBox(G2frame.dataWindow,choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        numclust.SetValue(str(ClusData['NumClust']))
        numclust.Bind(wx.EVT_COMBOBOX,OnClusNum)
        kmeanssizer.Add(numclust,0,WACV)
        compute = wx.Button(G2frame.dataWindow,label='Compute')
        compute.Bind(wx.EVT_BUTTON,OnCompute)
        kmeanssizer.Add(compute)
        return kmeanssizer

    def OnPlotSel(event):
        ClusData['plots'] = plotsel.GetValue()
        if ClusData['plots'] == 'Suprise':
            G2plt.PlotClusterXYZ(G2frame,None,None,ClusData,PlotName='Suprise')
        else:
            G2plt.PlotClusterXYZ(G2frame,YM,XYZ,ClusData,PlotName=ClusData['Method'],Title=ClusData['Method'])

    def ScikitSizer():

        def OnClusMethod(event):
            ClusData['Scikit'] = clusMethod.GetValue()
            OnCompute(event)

        def OnClusNum(event):
            ClusData['NumClust'] = int(numclust.GetValue())
            OnCompute(event)

        def OnCompute(event):
            whitMat = SCV.whiten(ClusData['DataMatrix'])
            if ClusData['Scikit'] == 'K-Means':
                result = SKC.KMeans(n_clusters=ClusData['NumClust'],algorithm='elkan',init='k-means++').fit(whitMat)
                print('K-Means sum squared dist. to means %.2f'%result.inertia_)
            elif ClusData['Scikit'] == 'Spectral clustering':
                result = SKC.SpectralClustering(n_clusters=ClusData['NumClust']).fit(whitMat)
            elif ClusData['Scikit'] == 'Mean-shift':
                result = SKC.MeanShift().fit(whitMat)
                print('Number of Mean-shift clusters found: %d'%(np.max(result.labels_)+1))
            elif ClusData['Scikit'] == 'Affinity propagation':
                result = SKC.AffinityPropagation(affinity='precomputed',damping=0.5).fit(SSD.squareform(ClusData['ConDistMat']))
                print('Number of Affinity propagation clusters found: %d'%(np.max(result.labels_)+1))
            elif ClusData['Scikit'] == 'Agglomerative clustering':
                result = SKC.AgglomerativeClustering(n_clusters=ClusData['NumClust'],
                    affinity='precomputed',linkage='average').fit(SSD.squareform(ClusData['ConDistMat']))

            ClusData['codes'] = result.labels_
            ClusData['Metrics'] = Metrics(whitMat,result)
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData)

        def Metrics(whitMat,result):
            if np.max(result.labels_) >= 1:
                Scoeff = SKM.silhouette_score(whitMat,result.labels_,metric='euclidean')
                print('Silhouette Coefficient: %.3f'%Scoeff)
                CHcoeff = SKM.calinski_harabasz_score(whitMat,result.labels_)
                print('Calinski-Harabasz index (Variance ratio): %.3f'%CHcoeff)
                DBcoeff = SKM.davies_bouldin_score(whitMat,result.labels_)
                print('Davies-Bouldin Index: %.3f'%DBcoeff)
                return Scoeff,CHcoeff,DBcoeff
            else:
                print('number of clusters found must be > 1 for metrics to be determined')
                return None

        scikitSizer = wx.BoxSizer(wx.VERTICAL)
        txt = wx.StaticText(G2frame.dataWindow,label=
                    'If you use Scikit-Learn Cluster Analysis, please cite: '+
                    G2G.GetCite('Scikit-Learn'))
        txt.Wrap(500)
        scikitSizer.Add(txt)
        choice = ['K-Means','Affinity propagation','Mean-shift','Spectral clustering','Agglomerative clustering']
        clusSizer = wx.BoxSizer(wx.HORIZONTAL)
        clusSizer.Add(wx.StaticText(G2frame.dataWindow,label='Select clustering method: '),0,WACV)
        clusMethod = wx.ComboBox(G2frame.dataWindow,choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        clusMethod.SetValue(ClusData['Scikit'])
        clusMethod.Bind(wx.EVT_COMBOBOX,OnClusMethod)
        clusSizer.Add(clusMethod,0,WACV)
        if ClusData['Scikit'] in ['K-Means','Spectral clustering','Agglomerative clustering']:
            nchoice = [str(i) for i in range(2,16)]
            clusSizer.Add(wx.StaticText(G2frame.dataWindow,label=' Select number of clusters (2-15): '),0,WACV)
            numclust = wx.ComboBox(G2frame.dataWindow,choices=nchoice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            numclust.SetValue(str(ClusData['NumClust']))
            numclust.Bind(wx.EVT_COMBOBOX,OnClusNum)
            clusSizer.Add(numclust,0,WACV)
        compute = wx.Button(G2frame.dataWindow,label='Compute')
        compute.Bind(wx.EVT_BUTTON,OnCompute)
        clusSizer.Add(compute,0,WACV)
        scikitSizer.Add(clusSizer)
        useTxt = '%s used the whitened data matrix'%ClusData['Scikit']
        if ClusData['Scikit'] in ['Agglomerative clustering','Affinity propagation']:
            useTxt = '%s used %s for distance method'%(ClusData['Scikit'],ClusData['Method'])
        scikitSizer.Add(wx.StaticText(G2frame.dataWindow,label=useTxt))
        if ClusData.get('Metrics',None) is not None:
            metrics = ClusData['Metrics']
            scikitSizer.Add(wx.StaticText(G2frame.dataWindow,
                label='Metrics: Silhoutte: %.3f, Variance: %.3f, Davies-Bouldin: %.3f'%(metrics[0],metrics[1],metrics[2])))
        return scikitSizer

    def memberSizer():

        def OnClusNum(event):
            shoNum = int(numclust.GetValue())
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData,shoNum)

        def OnSelection(event):
            name = cluslist.GetStringSelection()
            item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
            G2frame.PatternId = item
            if 'PWDR' in name:
                G2pwpl.PlotPatterns(G2frame,newPlot=False,plotType='PWDR')
            else: #PDF
                data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame, item,'PDF Controls'))
                G2plt.PlotISFG(G2frame,data,plotType='G(R)')

        #need 15 colors; values adjusted to match xkcs/PCA plot colors. NB: RGB reverse order from xkcd values.
        Colors = [['xkcd:blue',0xff0000],['xkcd:red',0x0000ff],['xkcd:green',0x00a000],['xkcd:cyan',0xd0d000],
                  ['xkcd:magenta',0xa000a0],['xkcd:black',0x000000],['xkcd:pink',0xb469ff],['xkcd:brown',0x13458b],
                  ['xkcd:teal',0x808000],['xkcd:orange',0x008cff],['xkcd:grey',0x808080],['xkcd:violet',0xe22b8a],
                  ['xkcd:aqua',0xaaaa00],['xkcd:blueberry',0xcd5a6a],['xkcd:bordeaux',0x00008b]]
        NClust = np.max(ClusData['codes'])
        memSizer = wx.BoxSizer(wx.VERTICAL)
        memSizer.Add(wx.StaticText(G2frame.dataWindow,label='Cluster populations (colors refer to cluster colors in PCA plot):'))
        for i in range(NClust+1):
            nPop= len(ClusData['codes'])-np.count_nonzero(ClusData['codes']-i)
            txt = wx.StaticText(G2frame.dataWindow,label='Cluster #%d has %d members'%(i,nPop))
            txt.SetForegroundColour(wx.Colour(Colors[i][1]))
            if wx.Colour(Colors[i][1]).GetLuminance() > 0.5:
                txt.SetBackgroundColour(wx.Colour(50,50,50))
            else:
                txt.SetBackgroundColour(wx.Colour(200,200,200))
            memSizer.Add(txt)
        headSizer = wx.BoxSizer(wx.HORIZONTAL)
        headSizer.Add(wx.StaticText(G2frame.dataWindow,label='Select cluster to list members: '),0,WACV)
        choice = [str(i) for i in range(NClust+1)]
        numclust = wx.ComboBox(G2frame.dataWindow,choices=choice,value=str(shoNum),style=wx.CB_READONLY|wx.CB_DROPDOWN)
        numclust.Bind(wx.EVT_COMBOBOX,OnClusNum)
        headSizer.Add(numclust,0,WACV)
        memSizer.Add(headSizer)
        if shoNum >= 0:
            memSizer.Add(wx.StaticText(G2frame.dataWindow,label='Members of cluster %d (select to show data plot):'%shoNum))
            ClusList = []
            for i,item in enumerate(ClusData['Files']):
                 if ClusData['codes'][i] == shoNum:
                     ClusList.append(item)
            cluslist = wx.ListBox(G2frame.dataWindow, choices=ClusList)
            cluslist.SetForegroundColour(wx.Colour(Colors[shoNum][1]))
            if wx.Colour(Colors[shoNum][1]).GetLuminance() > 0.5:
                cluslist.SetBackgroundColour(wx.Colour(50,50,50))
            else:
                cluslist.SetBackgroundColour(wx.Colour(200,200,200))
            cluslist.Bind(wx.EVT_LISTBOX,OnSelection)
            memSizer.Add(cluslist)
        return memSizer

    def outlierSizer():

        def OnOutSel(event):
            ClusData['OutMethod'] = outsel.GetValue()
            OnCompute(event)

        def OnCompute(event):
            if ClusData['OutMethod'] == 'One-Class SVM':
                ClusData['codes'] = SKVM.OneClassSVM().fit_predict(ClusData['DataMatrix'])  #codes = 1 or -1
            elif ClusData['OutMethod'] == 'Isolation Forest':
                ClusData['codes'] = SKE.IsolationForest().fit_predict(ClusData['DataMatrix'])
            elif ClusData['OutMethod'] == 'Local Outlier Factor':
                ClusData['codes'] = SKN.LocalOutlierFactor().fit_predict(ClusData['DataMatrix'])
            wx.CallAfter(UpdateClusterAnalysis,G2frame,ClusData,shoNum)

        outSizer = wx.BoxSizer(wx.VERTICAL)
        outSizer.Add(wx.StaticText(G2frame.dataWindow,label='Outlier (bad data) analysis with Scikit-learn:'))
        choice = ['One-Class SVM','Isolation Forest','Local Outlier Factor']
        outline = wx.BoxSizer(wx.HORIZONTAL)
        outline.Add(wx.StaticText(G2frame.dataWindow,label='Select method: '),0,WACV)
        outsel = wx.ComboBox(G2frame.dataWindow,choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
        outsel.SetValue(ClusData['OutMethod'])
        outsel.Bind(wx.EVT_COMBOBOX,OnOutSel)
        outline.Add(outsel,0,WACV)
        compute = wx.Button(G2frame.dataWindow,label='Compute')
        compute.Bind(wx.EVT_BUTTON,OnCompute)
        outline.Add(compute,0,WACV)
        outSizer.Add(outline)
        return outSizer

    def OnSelection(event):
        name = outlist.GetStringSelection()
        item = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
        G2frame.PatternId = item
        if 'PWDR' in name:
            G2frame.PatternId = G2gd.GetGPXtreeItemId(G2frame,G2frame.root,name)
            G2pwpl.PlotPatterns(G2frame,newPlot=False,plotType='PWDR')
        else:
            data = G2frame.GPXtree.GetItemPyData(G2gd.GetGPXtreeItemId(G2frame, item,'PDF Controls'))
            G2plt.PlotISFG(G2frame,data,plotType='G(R)')

    #patch
    ClusData['SKLearn'] = ClusData.get('SKLearn',False)
    ClusData['plots'] = ClusData.get('plots','All')
    ClusData['Scikit'] = ClusData.get('Scikit','K-Means')
    ClusData['OutMethod'] = ClusData.get('OutMethod','Isolation Forest')
    ClusData['MinkP'] = ClusData.get('MinkP','2')
    #end patch

    G2frame.dataWindow.ClearData()
    topSizer = G2frame.dataWindow.topBox
    topSizer.Clear(True)
    parent = G2frame.dataWindow.topPanel
    topSizer.Add(wx.StaticText(parent,label='Scipy Cluster Analysis:'),0,WACV)
    topSizer.Add((-1,-1),1,wx.EXPAND)
    topSizer.Add(G2G.HelpButton(parent,helpIndex=G2frame.dataWindow.helpKey))
    wx.CallAfter(G2frame.dataWindow.SetDataSize)
    bigSizer = wx.BoxSizer(wx.HORIZONTAL)
    mainSizer = wx.BoxSizer(wx.VERTICAL)
    mainSizer.Add((5,5),0)

    mainSizer.Add(FileSizer())
    if len(ClusData['Files']):
        mainSizer.Add(LimitSizer())
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,label='Cluster Analysis input data size: %d'%(GetYMatSize())))
        mainSizer.Add(wx.StaticText(G2frame.dataWindow,label=
            '(Examine any %s plot for reasonable limits; any change will clear Cluster data matrix) '%ClusData['Type']))
        makeArray = wx.Button(G2frame.dataWindow,label='Make Cluster Analysis data array')
        makeArray.Bind(wx.EVT_BUTTON,OnMakeArray)
        mainSizer.Add(makeArray)
        if len(ClusData['DataMatrix']):

            G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
            mainSizer.Add(wx.StaticText(G2frame.dataWindow,label='Distance Cluster Analysis:'))
            mainSizer.Add(MethodSizer())
            YM = None
            if len(ClusData['ConDistMat']):
                YM = SSD.squareform(ClusData['ConDistMat'])
                U,s,VT = nl.svd(YM) #s are the Eigenvalues
                ClusData['PCA'] = s
                s[3:] = 0.
                S = np.diag(s)
                XYZ = np.dot(S,VT)
                G2plt.PlotClusterXYZ(G2frame,YM,XYZ[:3,:],ClusData,PlotName=ClusData['Method'],Title=ClusData['Method'])
                G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
                mainSizer.Add(wx.StaticText(G2frame.dataWindow,label='Hierarchical Cluster Analysis:'))
                mainSizer.Add(HierSizer())

                G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
                mainSizer.Add(wx.StaticText(G2frame.dataWindow,label='K-means Cluster Analysis:'))
                mainSizer.Add(kmeanSizer())
                if 'dists' in ClusData:
                    kmeansres = wx.BoxSizer(wx.HORIZONTAL)
                    kmeansres.Add(wx.StaticText(G2frame.dataWindow,label='K-means ave. dist = %.2f'%np.mean(ClusData['dists'])))
                    mainSizer.Add(kmeansres)
            if ClusData['codes'] is not None:
                G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
                mainSizer.Add(memberSizer())
            G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
            plotSizer = wx.BoxSizer(wx.HORIZONTAL)
            plotSizer.Add(wx.StaticText(G2frame.dataWindow,label='Plot selection: '),0,WACV)
            if ClusData['CLuZ'] is None:
                choice = ['All','Distances','3D PCA','2D PCA','Diffs','Suprise']
            else:
                choice = ['All','Distances','Dendrogram','2D PCA','3D PCA','Diffs','Suprise']
            plotsel = wx.ComboBox(G2frame.dataWindow,choices=choice,style=wx.CB_READONLY|wx.CB_DROPDOWN)
            plotsel.SetValue(str(ClusData['plots']))
            plotsel.Bind(wx.EVT_COMBOBOX,OnPlotSel)
            plotSizer.Add(plotsel,0,WACV)
            mainSizer.Add(plotSizer)

            if ClusData['SKLearn'] and len(ClusData['ConDistMat']):
                G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
                subSizer = wx.BoxSizer(wx.HORIZONTAL)
                subSizer.Add((-1,-1),1,wx.EXPAND)
                subSizer.Add(wx.StaticText(G2frame.dataWindow,label='Scikit-Learn Cluster Analysis: '),0,WACV)
                subSizer.Add((-1,-1),1,wx.EXPAND)
                mainSizer.Add(subSizer,0,wx.EXPAND)
                mainSizer.Add(ScikitSizer())

        if ClusData['SKLearn'] and len(ClusData['DataMatrix']) > 15:
            G2G.HorizontalLine(mainSizer,G2frame.dataWindow)
            mainSizer.Add(outlierSizer())
            Nout = 0
            if ClusData['codes'] is not None:
                Nout = len(ClusData['codes'])-np.count_nonzero(ClusData['codes']+1)
            if Nout > 0:
                mainSizer.Add(wx.StaticText(G2frame.dataWindow,
                    label='%d Probable outlier data found by %s (select to show data plot):'%(Nout,ClusData['OutMethod'])))
                OutList = []
                for i,item in enumerate(ClusData['Files']):
                     if ClusData['codes'][i] < 0:
                         OutList.append(item)
                outlist = wx.ListBox(G2frame.dataWindow, choices=OutList)
                outlist.Bind(wx.EVT_LISTBOX,OnSelection)
                mainSizer.Add(outlist)
            else:
                mainSizer.Add(wx.StaticText(G2frame.dataWindow,label='No outlier data found'))

    bigSizer.Add(mainSizer)
    bigSizer.Layout()
    bigSizer.FitInside(G2frame.dataWindow)
    G2frame.dataWindow.SetSizer(bigSizer)
    G2frame.dataWindow.SetDataSize()
    G2frame.SendSizeEvent()
