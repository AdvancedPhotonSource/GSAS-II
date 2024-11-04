#testMagSym
import sys
import wx
import numpy as np
import copy
import GSASIIpath
GSASIIpath.SetBinaryPath()
import GSASIIspc as G2spc
import GSASIIctrlGUI as G2G
import GSASIIphsGUI as G2phsGUI

try:
    wx.NewId
except AttributeError:
    wx.NewId = wx.NewIdRef
    
[wxID_FILEEXIT, 
] = [wx.NewId() for _init_coll_File_Items in range(1)]
WACV = wx.ALIGN_CENTER_VERTICAL

class testMagSym(wx.Frame):

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='testMagSym', parent=parent,
            size=[800,300],style=wx.DEFAULT_FRAME_STYLE, title='Test magnetic space groups')
        self.testSSMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.File.Append(wxID_FILEEXIT,'Exit')
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        self.testSSMenu.Append(menu=self.File, title='Run')
        self.SetMenuBar(self.testSSMenu)
        self.testSSPanel = wx.Window(self)        
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)    
        self.dataFrame = None
        Data = {'SGData':G2spc.SpcGroup('P 1')[1],}
        self.propVec = '0,0,0'
        self.XYZ = '0,0,0'
        self.UpdateData(Data)

    def ExitMain(self, event):
        sys.exit()
        
    def OnFileExit(self,event):
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()
        
    def UpdateData(self,Data):
        
        def OnSpaceGroup(event):
            SpGrp = 'P 1'
            SpcGp = G2phsGUI.GetSpGrpfromUser(self,SpGrp)
            SGErr,SGData = G2spc.SpcGroup(SpcGp)
            if SGErr:
                text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                SGTxt.SetLabel(Data['SGData']['SpGrp'])
                msg = 'Space Group Error'
                Style = wx.ICON_EXCLAMATION
                Text = '\n'.join(text)
                wx.MessageBox(Text,caption=msg,style=Style)
            else:
                text,table = G2spc.SGPrint(SGData,AddInv=True)
                Data['SGData'] = SGData
                SGTxt.SetLabel(Data['SGData']['SpGrp'])
                msg = 'Space Group Information'
                SgNo = G2spc.SpaceGroupNumber(SpcGp)
                if SgNo+1:
                    text[0] += ', No. '+str(SgNo)
                G2G.SGMessageBox(self,msg,text,table).Show()
            self.UpdateData(Data)

        def OnSpinOp(event):
            Obj = event.GetEventObject()
            isym = Indx[Obj.GetId()]+1
            spCode = {'red':-1,'black':1}                    
            SGData['SGSpin'][isym] = spCode[Obj.GetValue()]
            G2spc.CheckSpin(isym,SGData)
            GenSym,GenFlg,BNSsym = G2spc.GetGenSym(SGData)
            SGData['GenSym'] = GenSym
            SGData['GenFlg'] = GenFlg
            SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
            wx.CallAfter(self.UpdateData,Data)
            
        def OnBNSlatt(event):
            Obj = event.GetEventObject()
            BNSlatt = Obj.GetValue()
            SGData.update(G2spc.SpcGroup(SGData['SpGrp'])[1])
#            SGData['SGSpin'] = [1,]*len(SGData['SGSpin'])
            if not BNSlatt:
                SGData['BNSlattsym'] = ['',[]]
            else:
                SGData['BNSlattsym'] = [BNSlatt,BNSsym[BNSlatt]]
            self.Trans = G2spc.ApplyBNSlatt(SGData,SGData['BNSlattsym'])
            GenSym,GenFlg = G2spc.GetGenSym(SGData)[:2]
            SGData['GenSym'] = GenSym
            SGData['GenFlg'] = GenFlg
            SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
            wx.CallAfter(self.UpdateData,Data)
            
        def OnXYZ(event):
            Obj = event.GetEventObject()
            self.XYZ = Obj.GetValue()
            try:
                XYZ= np.array(eval('['+self.XYZ+']'),dtype=np.float64)
                print('for:',XYZ)
                SytSym,Mul,Nop,dupDir = G2spc.SytSym(XYZ,SGData)
                print('dupDir',dupDir.keys())
                MagSytSym = G2spc.MagSytSym(SytSym,dupDir,SGData)
                CSI = G2spc.GetCSpqinel(SGData['SpnFlp'],dupDir)
#                print('for site sym = %s, mag site sym %s'%(SytSym.strip(),MagSytSym))
#                for opr in dupDir:
#                    print('opr: %s, no. %d, spin %d'%(opr,dupDir[opr],SGData['SpnFlp'][dupDir[opr]]))
                magStr = [str(smag) for smag in CSI[0]]
                magVal = [str(val) for val in CSI[1]]
                CSItxt.SetLabel(' mag. site sym: %6s, mult: %3d, CSI: %s %s'%(MagSytSym,Mul,magStr,magVal))
            except:
                print('Bad X,Y,Z entry: ',Obj.GetValue())
                self.XYZ = '0,0,0'
                Obj.SetValue(self.XYZ)
                
        def OnPropVec(event):
            Obj = event.GetEventObject()
            self.propVec = Obj.GetValue()
            try:
                pXYZ = np.array(eval('['+self.propVec+']'),dtype=np.float64)
                print('Little group operators for propagation vector:',pXYZ)
                Little = G2spc.GetLittleGrpOps(SGData,pXYZ)
                Lirt = G2spc.PackRot(Little)
                OpNames = [G2spc.GetOprName(str(irt)) for irt in Lirt]
                for Opr,name in zip(Little,OpNames):
                    print('  %s:  %s'%(G2spc.MT2text(Opr),name))
            except:
                print('Bad prop vector entry: ',Obj.GetValue())
                self.propVec = '0,0,0'
                Obj.SetValue(self.propVec)
                
        def OnShowOps(event):
            text,table = G2spc.SGPrint(SGData)
            msg = 'Space Group Information'
            SgNo = G2spc.SpaceGroupNumber(SGData['SpGrp'])
            if SgNo+1:
                text[0] += ', No. '+str(SgNo)
            G2G.SGMessageBox(self,msg,text,table,SpnFlp).Show()
            
        def OnShowGen(event):
            GenText = G2spc.TextGen(SGData,reverse=True)
            print(' Symmetry generators for %s:'%text[0].split(':')[1])
            for item in GenText:
                print(item)
                           
        def OnShowMOps(event):
            text,table = G2spc.SGPrint(SGData,AddInv=True)
            text[0] = ' Magnetic Space Group: '+SGData['MagSpGrp']
            text[3] = ' The magnetic lattice point group is '+SGData['MagPtGp']
            if SGData['SGGray'] and "1'" not in text[0]:
                text[0] += " 1'"
                text[3] += "1'"
            G2G.SGMagSpinBox(self.testSSPanel,msg,text,table,SGData['SGCen'],OprNames,SpnFlp,SGData['SGGray']).Show()

        def OnTestHKL(event):
            print('Extinctions for '+Data['SGData']['MagSpGrp'])
            hkls = np.mgrid[-6:6,-6:6,-6:6]
            HKLs = hkls.reshape((3,-1)).T
            for hkl in HKLs[1:]:    #skip 0,0,0
                ext = G2spc.checkMagextc(hkl,SGData)
                if ext: print(hkl)

        SGData = Data['SGData']
        Nops = len(SGData['SGOps'])*len(SGData['SGCen'])
        GenSym,GenFlg,BNSsym = G2spc.GetGenSym(SGData)
        Data['SGData']['GenSym'] = GenSym
        Data['SGData']['GenFlg'] = GenFlg
        if SGData['SGGray']:
            Data['SGData']['SpnFlp'] = Nops*[1,]+Nops*[-1,]
        SGData = Data['SGData']
        Data['SGData']['MagSpGrp'] = G2spc.MagSGSym(SGData)
        self.testSSPanel.DestroyChildren()
        mainSizer = wx.FlexGridSizer(0,2,5,5)
        mainSizer.Add(wx.StaticText(self.testSSPanel,-1,'  Space group: '),0,WACV)
        SpGrp = Data['SGData']['SpGrp']
        if SGData['SGGray']: SpGrp += " 1'"
        SGTxt = wx.Button(self.testSSPanel,wx.ID_ANY,SpGrp,size=(100,-1))
        SGTxt.Bind(wx.EVT_BUTTON,OnSpaceGroup)
        mainSizer.Add(SGTxt,0,WACV)
        mainSizer.Add(wx.StaticText(self.testSSPanel,-1,'  Propagation vector: '),0,WACV)
        propVec = wx.TextCtrl(self.testSSPanel,value=self.propVec,style=wx.TE_PROCESS_ENTER)
        propVec.Bind(wx.EVT_TEXT_ENTER,OnPropVec)        
        mainSizer.Add(propVec,0,WACV)
        magSizer = wx.BoxSizer(wx.VERTICAL)
        magSizer.Add(wx.StaticText(self.testSSPanel,label=' Magnetic spin operator selection:'))
        mainSizer.Add(magSizer)    
        spinSizer = wx.BoxSizer(wx.HORIZONTAL)
        Indx = {}
        spinColor = ['black','red']
        spCode = {-1:'red',1:'black'}
        for isym,sym in enumerate(SGData['GenSym'][1:]):
            spinSizer.Add(wx.StaticText(self.testSSPanel,label=' %s: '%(sym.strip())),0,WACV)                
            spinOp = wx.ComboBox(self.testSSPanel,value=spCode[SGData['SGSpin'][isym+1]],choices=spinColor,
                style=wx.CB_READONLY|wx.CB_DROPDOWN)                
            Indx[spinOp.GetId()] = isym
            spinOp.Bind(wx.EVT_COMBOBOX,OnSpinOp)
            spinSizer.Add(spinOp,0,WACV)
        if '(' in SGData['BNSlattsym'][0]:
            spinSizer.Add(wx.StaticText(self.testSSPanel,label='Choose new space group to change BNS'),0,WACV)
        else:
            if BNSsym:
                spinSizer.Add(wx.StaticText(self.testSSPanel,label=' BNS lattice:'),0,WACV)
                BNS = wx.ComboBox(self.testSSPanel,value=SGData['BNSlattsym'][0],choices=['',] +list(BNSsym.keys()),style=wx.CB_READONLY|wx.CB_DROPDOWN)
                BNS.Bind(wx.EVT_COMBOBOX,OnBNSlatt)
                spinSizer.Add(BNS,0,WACV)
        mainSizer.Add(spinSizer)
        OprNames,SpnFlp = G2spc.GenMagOps(SGData)
        SGData['SpnFlp'] = SpnFlp
        SGData['OprNames'] = OprNames
        msg = 'Magnetic space group information'
        text,table = G2spc.SGPrint(SGData,AddInv=True)
        text[0] = ' Magnetic Space Group: '+SGData['MagSpGrp']
        text[3] = ' The magnetic lattice point group is '+SGData['MagPtGp']
        if SGData['SGGray'] and "1'" not in text[0]:
            text[0] += " 1'"
            text[3] += "1'"
        mainSizer.Add(wx.StaticText(self.testSSPanel,label=' Magnetic space group: %s  '%(SGData['MagSpGrp'])),0,WACV)
        mainSizer.Add(wx.StaticText(self.testSSPanel,label='Mag Gen: %s'%str(SGData['SGSpin'])))
        SpPos = wx.BoxSizer(wx.HORIZONTAL)
        SpPos.Add(wx.StaticText(self.testSSPanel,label=' X,Y,Z:'),0,WACV)
        xyz = wx.TextCtrl(self.testSSPanel,value=self.XYZ,style=wx.TE_PROCESS_ENTER)
        xyz.Bind(wx.EVT_TEXT_ENTER,OnXYZ)        
        SpPos.Add(xyz,0,WACV)
        mainSizer.Add(SpPos,0,WACV)
        XYZ= np.array(eval('['+self.XYZ+']'),dtype=np.float64)
        SytSym,Mul,Nop,dupDir = G2spc.SytSym(XYZ,SGData)
        CSI = G2spc.GetCSpqinel(SGData['SpnFlp'],dupDir)
        magStr = [str(smag) for smag in CSI[0]]
        magVal = [str(val) for val in CSI[1]]
        CSItxt = wx.StaticText(self.testSSPanel,label=' site sym: %6s, mult: %3d, MCSI: %s %s'%(SytSym,Mul,magStr,magVal))
        mainSizer.Add(CSItxt,0,WACV)
        testHKL = wx.Button(self.testSSPanel,-1,'Extinction test')
        testHKL.Bind(wx.EVT_BUTTON,OnTestHKL)
        mainSizer.Add(testHKL,0,WACV)
        printSizer = wx.BoxSizer(wx.HORIZONTAL)
        showOps = wx.Button(self.testSSPanel,-1,'Show sym. ops')
        showOps.Bind(wx.EVT_BUTTON,OnShowOps)
        printSizer.Add(showOps,0,WACV)
        showGen = wx.Button(self.testSSPanel,-1,'Print generators')
        showGen.Bind(wx.EVT_BUTTON,OnShowGen)
        printSizer.Add(showGen,0,WACV)
        showMOps = wx.Button(self.testSSPanel,-1,'Show mag. ops')
        showMOps.Bind(wx.EVT_BUTTON,OnShowMOps)
        printSizer.Add(showMOps,0,WACV)
        mainSizer.Add(printSizer,0,WACV)
        SGData1 = copy.deepcopy(SGData)
        SGData1['SGSpin'] = G2spc.GetSGSpin(SGData1,SGData1['MagSpGrp'])
        mainSizer.Add(wx.StaticText(self.testSSPanel,label='New Mag Gen: %s'%str(SGData1['SGSpin'])))
        SGData1['GenSym'],SGData1['GenFlg'],BNSsym = G2spc.GetGenSym(SGData1)
        MagSpGrp = G2spc.MagSGSym(SGData1).replace(' ','')
        mainSizer.Add(wx.StaticText(self.testSSPanel,label='Gives symbol: %s'%MagSpGrp))
        
        SpGp = MagSpGrp.replace("'",'')
        SpGrp = G2spc.StandardizeSpcName(SpGp)
        SGData1 = G2spc.SpcGroup(SpGrp)[1]
        MSpGrp = G2spc.SplitMagSpSG(MagSpGrp)
        SGData1['SGSpin'] = G2spc.GetSGSpin(SGData1,MSpGrp)
        SGData1['GenSym'],SGData1['GenFlg'],BNSsym = G2spc.GetGenSym(SGData1)
        if '_' in MSpGrp[0]:
            SGData1['BNSlattsym'] = [MSpGrp[0],BNSsym[MSpGrp[0]]]
        else:
            SGData1['BNSlattsym'] = [SGData['SGLatt'],[0.,0.,0.]]
        mainSizer.Add(wx.StaticText(self.testSSPanel,label='Symbol Mag Gen: %s'%str(SGData1['SGSpin'])))
        MagSpGrp = G2spc.MagSGSym(SGData1).replace(' ','')
        mainSizer.Add(wx.StaticText(self.testSSPanel,label='Gives symbol: %s'%MagSpGrp))
        self.testSSPanel.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.testSSPanel)
        Size[0] = 800
        Size[1] = max(Size[1],350)
        self.testSSPanel.SetSize(Size)
        G2G.SGMagSpinBox(self.testSSPanel,msg,text,table,SGData['SGCen'],OprNames,SpnFlp,SGData['SGGray']).Show()
            
class testMagSmain(wx.App):
    def OnInit(self):
        self.main = testMagSym(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True

def main():
    'Starts main application to compute and plot derivatives'
#    GSASIIpath.InvokeDebugOpts()
    application = testMagSmain(0)
    application.MainLoop()
    
if __name__ == '__main__':
    main()
            
        
