#test
import sys
import os
import wx
import importlib.util
import numpy as np
try:
    importlib.util.find_spec('GSASII.GSASIIGUI')
except ModuleNotFoundError:
    print('GSAS-II not installed in Python; Hacking sys.path')
    sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
from GSASII import GSASIIpath
GSASIIpath.SetBinaryPath()
from GSASII import GSASIIspc as G2spc
from GSASII import GSASIIlattice as G2lat
from GSASII import GSASIIctrlGUI as G2G

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

[wxID_FILEEXIT, 
] = [wx.NewId() for _init_coll_File_Items in range(1)]
WACV = wx.ALIGN_CENTER_VERTICAL
laueSS = {'2/m':['(a0g)','(a1/2g)','(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)'],
        'mmm':['(a00)','(a1/20)','(a01/2)','(a1/21/2)','(a10)','(a01)',
               '(0b0)','(1/2b0)','(0b1/2)','(1/2b1/2)','(1b0)','(0b1)',
               '(00g)','(1/20g)','(01/2g)','(1/21/2g)','(10g)','(01g)']}
        
laueTS = {'2/m':['','s','s0','0s','ss'],
          'mmm':['000','s00','0s0','00s','ss0','s0s','0ss','q00','0q0','00q','0qq','q0q','qq0'],
          }

class testSSymbols(wx.Frame):

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='testSSymbols', parent=parent,
            size=wx.Size(500,400),style=wx.DEFAULT_FRAME_STYLE, title='Test SS symbols')
        self.testSSMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.File.Append(wxID_FILEEXIT,'Exit','Exit from testSS')
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        self.testSSMenu.Append(menu=self.File, title='Run')
        self.SetMenuBar(self.testSSMenu)
        self.testSSPanel = wx.Window(self)        
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)    
        self.dataFrame = None
        self.XYZ = '0,0,0'
        self.waveOrder = 1
        Data = {'SGData':G2spc.SpcGroup('P 1')[1],'SuperSg':'(abg)',}
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
            Flds = SGTxt.GetValue().split()
            #get rid of extra spaces between fields first
            for fld in Flds: fld = fld.strip()
            SpcGp = ' '.join(Flds)
            # try a lookup on the user-supplied name
            SpGrpNorm = G2spc.StandardizeSpcName(SpcGp)
            if SpGrpNorm:
                SGErr,SGData = G2spc.SpcGroup(SpGrpNorm)
            else:
                SGErr,SGData = G2spc.SpcGroup(SpcGp)
            if SGErr:
                text = [G2spc.SGErrors(SGErr)+'\nSpace Group set to previous']
                SGTxt.SetValue(Data['SGData']['SpGrp'])
                msg = 'Space Group Error'
                Style = wx.ICON_EXCLAMATION
                Text = '\n'.join(text)
                wx.MessageBox(Text,caption=msg,style=Style)
            else:
                text,table = G2spc.SGPrint(SGData)
                Data['SGData'] = SGData
                SGTxt.SetValue(Data['SGData']['SpGrp'])
                msg = 'Space Group Information'
                G2G.SGMessageBox(self,msg,text,table).Show()
            SSChoice = G2spc.SSChoice(Data['SGData'])
            Data['SuperSg'] = SSChoice[0]
            self.UpdateData(Data)

        def OnSuperGp(event):
            SSymbol = superGp.GetValue()
            SpGrp = Data['SGData']['SpGrp']
            if Data['SGData']['SGGray']: SpGrp += " 1'"
            print('Try: %s%s'%(Data['SGData']['SpGrp'],SSymbol))
            E,SSGData = G2spc.SSpcGroup(Data['SGData'],SSymbol)
            if SSGData:
                text,table = G2spc.SSGPrint(Data['SGData'],SSGData)
                Data['SSGData'] = SSGData
                Data['SuperSg'] = SSymbol
                msg = 'Superspace Group Information'
                G2G.SGMessageBox(self,msg,text,table).Show()
                print('Super spacegroup operators for '+SSGData['SSpGrp'])
                for Op in SSGData['SSGOps']:
                    print(G2spc.SSMT2text(Op).replace(' ',''))
                if SGData['SGInv']:                                 
                    for Op in SSGData['SSGOps']:
                        Op = [-Op[0],-Op[1]%1.]
                        print(G2spc.SSMT2text(Op).replace(' ',''))
            else:
                text = [E+'\nSuperspace Group set to previous']
                superGp.SetValue(Data['SuperSg'])
                msg = 'Superspace Group Error'
                Style = wx.ICON_EXCLAMATION
                Text = '\n'.join(text)
                wx.MessageBox(Text,caption=msg,style=Style)
            wx.CallAfter(self.UpdateData,Data)
        
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
                XYZ= np.array(eval('['+self.XYZ+']'),dtype=float)
                print('for:',XYZ)

            except:
                print('Bad X,Y,Z entry: ',Obj.GetValue())
                self.XYZ = '0,0,0'
                Obj.SetValue(self.XYZ)
            wx.CallAfter(self.UpdateData,Data)
                
        def siteSymSizer(XYZ,SGData):
            SSsizer = wx.BoxSizer(wx.VERTICAL)
            SpPos = wx.BoxSizer(wx.HORIZONTAL)
            SpPos.Add(wx.StaticText(self.testSSPanel,label=' Special position test:'),0,WACV)
            SpPos.Add(wx.StaticText(self.testSSPanel,label=' X,Y,Z:'),0,WACV)
            xyz = wx.TextCtrl(self.testSSPanel,value=self.XYZ,style=wx.TE_PROCESS_ENTER)
            xyz.Bind(wx.EVT_TEXT_ENTER,OnXYZ)        
            SpPos.Add(xyz,0,WACV)
            XYZ= np.array(eval('['+self.XYZ+']'),dtype=float)
            SytSym,Mul,Nop,dupDir = G2spc.SytSym(XYZ,SGData)
            SpPos.Add(wx.StaticText(self.testSSPanel,label=' site sym: %6s, mult: %3d'%(SytSym,Mul)),0,WACV)
            SSsizer.Add(SpPos)
            CSIX = G2spc.GetCSxinel(SytSym)
            CSIU = G2spc.GetCSuinel(SytSym)
            CSIM = G2spc.GetCSpqinel(SGData['SpnFlp'],dupDir)
            StrXYZ = [str(sxyz) for sxyz in CSIX[0]]
            ValXYZ = [str(val) for val in CSIX[1]]
            StrUIJ = [str(suij) for suij in CSIU[0]]
            ValUIJ = [str(val) for val in CSIU[1]]
            StrM = [str(smag) for smag in CSIM[0]]
            ValM = [str(val) for val in CSIM[1]]
            CSIXtxt = wx.StaticText(self.testSSPanel,label=' CSI-X: %s %s'%(StrXYZ,ValXYZ))
            SSsizer.Add(CSIXtxt)
            SSsizer.Add((5,5),0)
            CSIUtxt = wx.StaticText(self.testSSPanel,label=' CSI-U: %s %s'%(StrUIJ,ValUIJ))
            SSsizer.Add(CSIUtxt)
            SSsizer.Add((5,5),0)
            CSIMtxt = wx.StaticText(self.testSSPanel,label=' CSI-M: %s %s'%(StrM,ValM))
            SSsizer.Add(CSIMtxt)
            SSsizer.Add((5,5),0)
            return SSsizer
        
        def waveSymSizer(XYZ,SGData,SSGData):
            
            def OnWave(event):
                self.waveOrder = int(wchoice.GetValue())
                wx.CallAfter(self.UpdateData,Data)
                
            sctxt = ['sin','cos']
            SSsizer = wx.BoxSizer(wx.VERTICAL)
            SStop = wx.BoxSizer(wx.HORIZONTAL)
            SStop.Add(wx.StaticText(self.testSSPanel,
                label=' Special position fourier wave %d test: Select wave '%self.waveOrder),0,WACV)
            wchoice = wx.ComboBox(self.testSSPanel,value=str(self.waveOrder),choices=['1','2','3','4'],style=wx.CB_DROPDOWN)
            wchoice.Bind(wx.EVT_COMBOBOX,OnWave)
            SStop.Add(wchoice,0,WACV)
            SSsizer.Add(SStop)
            for Stype in ['Spos','Sadp','Smag']:
                if self.waveOrder > 1 and Stype == 'Smag':
                    continue
                CSIW = G2spc.GetSSfxuinel('Fourier',Stype,self.waveOrder,XYZ,SGData,SSGData,True)[0]
                for i in [0,1]:
                    num = len(CSIW[0])//2
                    iBeg = i*num
                    iFin = (i+1)*num
                    StrSS = [str(ss[0]) for ss in CSIW[0][iBeg:iFin]]
                    ValSS = [str(val[0]) for val in CSIW[1][iBeg:iFin]]
                    CSIWtxt = wx.StaticText(self.testSSPanel,label=' CSI-%s %s: %s %s'%(Stype,sctxt[i],StrSS,ValSS))
                    SSsizer.Add(CSIWtxt)
                    SSsizer.Add((5,5),0)                
            return SSsizer
            

        SGData = Data['SGData']
        Nops = len(SGData['SGOps'])*len(SGData['SGCen'])
        GenSym,GenFlg,BNSsym = G2spc.GetGenSym(SGData)
        Data['SGData']['GenSym'] = GenSym
        Data['SGData']['GenFlg'] = GenFlg
        if SGData['SGGray']:
            Data['SGData']['SpnFlp'] = Nops*[1,]+Nops*[-1,]
        
        SSGData = G2spc.SSpcGroup(SGData,Data['SuperSg'])[1]
        self.testSSPanel.DestroyChildren()
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer = wx.FlexGridSizer(0,2,5,5)
        topSizer.Add(wx.StaticText(self.testSSPanel,-1,'  Space group: '),0,WACV)
        SpGrp = Data['SGData']['SpGrp']
        MagSpGrp = SGData['MagSpGrp'] = G2spc.MagSGSym(SGData)
        if Data['SGData']['SGGray']: SpGrp += " 1'"
        SGTxt = wx.TextCtrl(self.testSSPanel,-1,value=SpGrp,style=wx.TE_PROCESS_ENTER)
        SGTxt.Bind(wx.EVT_TEXT_ENTER,OnSpaceGroup)
        topSizer.Add(SGTxt,0,WACV)
        topSizer.Add(wx.StaticText(self.testSSPanel,label=' Superspace group: '+MagSpGrp),0,WACV)
        ssChoice = G2spc.SSChoice(Data['SGData'])
        ssSym = Data['SuperSg']
        if ssChoice:
            superGp = wx.ComboBox(self.testSSPanel,value=ssSym,choices=ssChoice,style=wx.CB_DROPDOWN)   #wx.CB_READONLY|
            superGp.Bind(wx.EVT_COMBOBOX,OnSuperGp)
            superGp.Bind(wx.EVT_TEXT_ENTER,OnSuperGp)
        else:   #nonstandard space group symbol not in my dictionary
            superGp = wx.TextCtrl(self.testSSPanel,value=ssSym,style=wx.TE_PROCESS_ENTER)
            superGp.Bind(wx.EVT_TEXT_ENTER,OnSuperGp)                        
        topSizer.Add(superGp,0,WACV)
        mainSizer.Add(topSizer)
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
        XYZ= np.array(eval('['+self.XYZ+']'),dtype=float)
        mainSizer.Add(siteSymSizer(XYZ,SGData))
        mainSizer.Add(waveSymSizer(XYZ,SGData,SSGData))
        
        self.testSSPanel.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.testSSPanel)
        Size[0] = max(Size[0],500)
        Size[1] = max(Size[1],400)
        self.testSSPanel.SetSize(Size)
        msg = 'Magnetic space group information'
        text,table = G2spc.SGPrint(SGData,AddInv=True)
        text[0] = ' Magnetic Space Group: '+SGData['MagSpGrp']
        text[3] = ' The magnetic lattice point group is '+SGData['MagPtGp']
        if SGData['SGGray'] and "1'" not in text[0]:
            text[0] += " 1'"
            text[3] += "1'"
        G2G.SGMagSpinBox(self.testSSPanel,msg,text,table,SGData['SGCen'],OprNames,SpnFlp,SGData['SGGray']).Show()
            
class testSSmain(wx.App):
    def OnInit(self):
        self.main = testSSymbols(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True

def main():
    'Starts main application to compute and plot derivatives'
    application = testSSmain(0)
    application.MainLoop()
    
if __name__ == '__main__':
    main()
            
        
