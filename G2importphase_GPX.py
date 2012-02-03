# Routines to import Phase information from GSAS-II .gpx files
import cPickle
import GSASIIIO as G2IO

class PhaseReaderClass(G2IO.ImportPhase):
    def __init__(self):
        super(self.__class__,self).__init__( # fancy way to say ImportPhase.__init__
            extensionlist=('.gpx',),
            strictExtension=True,
            formatName = 'GSAS-II gpx',
            longFormatName = 'GSAS-II project (.gpx file) import'
            )
    def ContentsValidator(self, filepointer):
        filepointer.seek(0) # rewind the file pointer
        # if the 1st section can't be read as a cPickle file, it can't be!
        try: 
            cPickle.load(filepointer)
        except:
            return False
        return True
    def Reader(self,filename,filepointer, ParentFrame=None):
        try:
            phasenames = G2IO.GetPhaseNames(filename)
            print phasenames
        except:
            return False
        if not phasenames:
            return False            # no blocks with coordinates
        elif len(phasenames) == 1: # no choices
            selblk = 0
        else:                       # choose from options                
            selblk = self.PhaseSelector(
                phasenames,
                ParentFrame=ParentFrame,
                title= 'Select a phase from the list below',
                )
            if selblk is None: return False # User pressed cancel
        try:
            self.Phase = G2IO.GetAllPhaseData(filename,phasenames[selblk])
            return True
        except:
            return False

