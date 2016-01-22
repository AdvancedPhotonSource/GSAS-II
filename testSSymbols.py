#test
import sys
import wx
import GSASIIspc as G2spc
import GSASIIgrid as G2gd
import numpy as np
import copy

def create(parent):
    return testDeriv(parent)
    
[wxID_FILEEXIT, 
] = [wx.NewId() for _init_coll_File_Items in range(1)]
WACV = wx.ALIGN_CENTER_VERTICAL

class testSSymbols(wx.Frame):

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='testSSymbols', parent=parent,
            size=wx.Size(300, 150),style=wx.DEFAULT_FRAME_STYLE, title='Test SS symbols')
        self.testSSMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.File.Append(help='Exit from testSS', id=wxID_FILEEXIT, kind=wx.ITEM_NORMAL,
            text='Exit')
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        self.testSSMenu.Append(menu=self.File, title='Run')
        self.SetMenuBar(self.testSSMenu)
        self.testSSPanel = wx.Window(self)        
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)    
        self.dataFrame = None
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
        
        def OnTryAll(event):
            SSList = G2spc.ssdict.get(Data['SGData']['SpGrp'],['',])
            for SSymbol in SSList:
                E,SSGData = G2spc.SSpcGroup(Data['SGData'],SSymbol)
                if SSGData:
                    text,table = G2spc.SSGPrint(Data['SGData'],SSGData)
                    Data['SSGData'] = SSGData
                    Data['SuperSg'] = SSymbol
                    msg = 'Superspace Group Information'
                    G2gd.SGMessageBox(self,msg,text,table).Show()
                else:
                    msg = 'Superspace Group Error for'+SSymbol
                    Style = wx.ICON_EXCLAMATION
                    Text = '\n'+E
                    wx.MessageBox(Text,caption=msg,style=Style)
        
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
                G2gd.SGMessageBox(self,msg,text,table).Show()
            SSChoice = G2spc.ssdict.get(Data['SGData']['SpGrp'],['',])
            Data['SuperSg'] = SSChoice[0]
            self.UpdateData(Data)

        def OnSuperGp(event):
            SSymbol = superGp.GetValue()
            E,SSGData = G2spc.SSpcGroup(Data['SGData'],SSymbol)
            if SSGData:
                text,table = G2spc.SSGPrint(Data['SGData'],SSGData)
                Data['SSGData'] = SSGData
                Data['SuperSg'] = SSymbol
                msg = 'Superspace Group Information'
                G2gd.SGMessageBox(self,msg,text,table).Show()
                print 'Super spacegroup operators for '+SSGData['SSpGrp']
                for Op in SSGData['SSGOps']:
                    print G2spc.SSMT2text(Op).replace(' ','')
                if SGData['SGInv']:                                 
                    for Op in SSGData['SSGOps']:
                        Op = [-Op[0],-Op[1]%1.]
                        print G2spc.SSMT2text(Op).replace(' ','')                                 
            else:
                text = [E+'\nSuperspace Group set to previous']
                superGp.SetValue(Data['SuperSg'])
                msg = 'Superspace Group Error'
                Style = wx.ICON_EXCLAMATION
                Text = '\n'.join(text)
                wx.MessageBox(Text,caption=msg,style=Style)
            self.UpdateData(Data)
        
        SGData = G2spc.SpcGroup(Data['SGData']['SpGrp'])[1]
        SSGData = G2spc.SSpcGroup(SGData,Data['SuperSg'])[1]
        
        self.testSSPanel.DestroyChildren()
        mainSizer = wx.FlexGridSizer(0,2,5,5)
        mainSizer.Add(wx.StaticText(self.testSSPanel,-1,'  Space group: '),0,WACV)
        SGTxt = wx.TextCtrl(self.testSSPanel,-1,value=Data['SGData']['SpGrp'],style=wx.TE_PROCESS_ENTER)
        SGTxt.Bind(wx.EVT_TEXT_ENTER,OnSpaceGroup)
        mainSizer.Add(SGTxt,0,WACV)
        mainSizer.Add(wx.StaticText(self.testSSPanel,label=' Superspace group: '+Data['SGData']['SpGrp']),0,WACV)
        SSChoice = G2spc.ssdict.get(Data['SGData']['SpGrp'],[])
        if SSChoice:
            superGp = wx.ComboBox(self.testSSPanel,value=Data['SuperSg'],choices=SSChoice,style=wx.CB_DROPDOWN)   #wx.CB_READONLY|
            superGp.Bind(wx.EVT_COMBOBOX,OnSuperGp)
            superGp.Bind(wx.EVT_TEXT_ENTER,OnSuperGp)
        else:   #nonstandard space group symbol not in my dictionary
            superGp = wx.TextCtrl(self.testSSPanel,value=Data['SuperSg'],style=wx.TE_PROCESS_ENTER)
            superGp.Bind(wx.EVT_TEXT_ENTER,OnSuperGp)                        
        mainSizer.Add(superGp,0,WACV)
        mainSizer.Add(wx.StaticText(self.testSSPanel,-1,' Try all SS symbols: '),0,WACV)
        SStry = wx.Button(self.testSSPanel,-1,'OK') 
        SStry.Bind(wx.EVT_BUTTON,OnTryAll)
        mainSizer.Add(SStry,0,WACV)
        self.testSSPanel.SetSizer(mainSizer)
        Size = mainSizer.Fit(self.testSSPanel)
        Size[0] = 800
        Size[1] = max(Size[1],350)
        self.testSSPanel.SetSize(Size)
            
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
            
        
