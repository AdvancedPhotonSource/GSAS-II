# -*- coding: utf-8 -*-
'''Classes in :mod:`~GSASII.exports.G2export_PDB` follow:
'''
from __future__ import division, print_function
import numpy as np
import os.path
from .. import GSASIIfiles as G2fil
from .. import GSASIIlattice as G2lat

class ExportPhasePDB(G2fil.ExportBaseclass):
    '''Used to create a PDB file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'PDB',
            extension='.PDB',
            longFormatName = 'Export phase as .PDB file'
            )
        self.exporttype = ['phase']
        self.multiple = True

    def Exporter(self,event=None):
        '''Export as a PDB file
        '''
        
        def PDBheader():
            self.Write("HEADER phase "+str(phasenam)+" from "+str(self.G2frame.GSASprojectfile))
            self.Write("TITLE")
            self.Write("COMPND")
            self.Write("SOURCE")
            self.Write("KEYWDS")
            self.Write("EXPDTA    X-RAY POWDER DIFFRACTION")
            self.Write("REVDAT")
            self.Write("JRNL")
            self.Write("REMARK   1")
            self.Write("REMARK   2")                                                                      
            self.Write("REMARK   2 RESOLUTION. 2.66 ANGSTROMS.")                                          
            self.Write("REMARK   2 POWDER DIFFRACTION MINIMUM D-SPACING.")
            
        def PDBremark250():                                
            self.Write('REMARK 250')                                                                      
            self.Write('REMARK 250 REFINEMENT.')                                                          
            self.Write('REMARK 250   PROGRAM     : GSAS-II')                                                 
            self.Write('REMARK 250   AUTHORS     : TOBY & VON DREELE')
            self.Write('REMARK 250   REFRENCE    : J. APPL. CRYST. 46, 544-549(2013)')                              
            self.Write('REMARK 250')
            self.Write('REMARK 250  DATA USED IN REFINEMENT')                                             
            self.Write('REMARK 250   RESOLUTION RANGE HIGH (ANGSTROMS) :  x.xx')                          
            self.Write('REMARK 250   RESOLUTION RANGE LOW  (ANGSTROMS) : xx.xx')                          
            self.Write('REMARK 250   POWDER DIFFRACTION DATA.')                                           
            self.Write('REMARK 250')
            self.Write('REMARK 250  FIT TO DATA USED IN REFINEMENT')                                      
            self.Write('REMARK 250   NUMBER OF POWDER PATTERNS         :     x')                          
            self.Write('REMARK 250   PROFILE R VALUES              (%) :  x.xx')                          
            self.Write('REMARK 250   WEIGHTED PROFILE R VALUES     (%) :  x.xx')                         
            self.Write('REMARK 250   F**2 R VALUES                 (%) : xx.xx')                          
            self.Write('REMARK 250   NUMBERS OF POWDER PATTERN POINTS  :  xxxx')                          
            self.Write('REMARK 250   NUMBERS OF REFLECTIONS            :  xxxx')                          
            self.Write('REMARK 250   TOTAL NUMBER OF POWDER POINTS     :  xxxx')                          
            self.Write('REMARK 250')
            self.Write('REMARK 250  NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.')                    
            self.Write('REMARK 250   PROTEIN ATOMS       :      xxxx')                                              
            self.Write('REMARK 250   NUCLEIC ACID ATOMS  :      xxxx')                                              
            self.Write('REMARK 250   HETEROGEN ATOMS     :      xxxx')                                              
            self.Write('REMARK 250   SOLVENT ATOMS       :      xxxx')                                             
            self.Write('REMARK 250')
            self.Write('REMARK 250  MODEL REFINEMENT.')                                                   
            self.Write('REMARK 250   NUMBER OF LEAST-SQUARES PARAMETERS :  xxxx')                         
            self.Write('REMARK 250   NUMBER OF RESTRAINTS               :  xxxx')                         
            self.Write('REMARK 250')
            self.Write('REMARK 250  RMS DEVIATIONS FROM RESTRAINT TARGET VALUES. NUMBER.')                
            self.Write('REMARK 250   BOND ANGLES                      (DEG) : x.xx   xxx')                
#            self.Write('REMARK 250   ANTI-BUMPING DISTANCE RESTRAINTS   (A) :x.xxx   xxx')                
#            self.Write('REMARK 250   HYDROGEN BOND DISTANCE RESTRAINTS  (A) :x.xxx   xxx')                
            self.Write('REMARK 250   INTERATOMIC DISTANCES              (A) :x.xxx   xxx')               
            self.Write('REMARK 250   DISTANCES FROM RESTRAINT PLANES    (A) :x.xxx   xxx')                
            self.Write('REMARK 250   TORSION PSEUDOPOTENTIAL RESTRAINTS (E) : x.xx   xxx')               
            self.Write('REMARK 250   TORSION ANGLE RESTRAINTS           (E) : x.xx   xxx')                
            self.Write('REMARK 250')
            self.Write('REMARK 200')                                                                      
            self.Write('DBREF')
            
        def PDBseqres(seqList):
            chains = list(seqList.keys())
            chains.sort()
            nSeq = 0
            for chain in chains:
                nres = len(seqList[chain])
                nrec = (nres-1)//13+1
                iB = 0
                for irec in range(nrec):
                    iF = min(iB+13,nres)
                    text = 'SEQRES {:3d}{:2s}{:5d}  '+(iF-iB)*'{:4s}'
                    self.Write(text.format(irec+1,chain,nres,*seqList[chain][iB:iF]))
                    nSeq += 1
                    iB += 13
            return nSeq
            
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()
        if self.ExportSelect():    # set export parameters; ask for file name
            return
        filename = self.filename
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info
            General = phasedict['General']
            if General['Type'] != 'macromolecular':
                print ('phase '+phasenam+' not macromolecular, skipping')
                continue
            i = self.Phases[phasenam]['pId']
            if len(self.phasenam) > 1: # if more than one filename is included, add a phase #
                self.filename = os.path.splitext(filename)[1] + "_" + str(i) + self.extension
            #fp = self.OpenFile()
            Atoms = phasedict['Atoms']
            cx,ct,cs,cia = General['AtomPtrs']
            seqList = {}
            AA3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
                'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE']
            seq = 0
            notProt = True
            for atom in Atoms:
                if atom[ct-3] in AA3letter and int(atom[ct-4]) != seq:
                    notProt = False
                    if atom[ct-2] not in seqList:
                        seqList[atom[ct-2]] = []
                    seqList[atom[ct-2]].append(atom[ct-3])
                    seq = int(atom[ct-4])
            PDBheader()
            PDBremark250()
            nSeq = PDBseqres(seqList)
            
            # get cell parameters
            Cell = General['Cell'][1:7]
            line = "CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} ".format(*Cell)
            line += General['SGData']['SpGrp'].ljust(13)
            line += '%2d'%(len(General['SGData']['SGOps'])*len(General['SGData']['SGCen']))
            self.Write(line)
            self.Write('ORIGX1      1.000000  0.000000  0.000000        0.00000')
            self.Write('ORIGX2      0.000000  1.000000  0.000000        0.00000')
            self.Write('ORIGX3      0.000000  0.000000  1.000000        0.00000')
            A,B = G2lat.cell2AB(Cell)
            self.Write('SCALE1     {:9.6f} {:9.6f} {:9.6f}        0.00000'.format(*B[0]))
            self.Write('SCALE2     {:9.6f} {:9.6f} {:9.6f}        0.00000'.format(*B[1]))
            self.Write('SCALE3     {:9.6f} {:9.6f} {:9.6f}        0.00000'.format(*B[2]))
            iatom = 1
            nHet = 0
            nTer = 0
            fmt = '{:6s}{:5d}  {:4s}{:3s} {:1s}{:4s}    '+3*'{:8.3f}'+2*'{:6.2f}'+'{:s}'
            for atom in Atoms:
                if atom[cia] == 'I':    #need to deal with aniso thermals for proteins = "ANISOU" records
                    Biso = atom[cia+1]*8.*np.pi**2
                xyz = np.inner(A,np.array(atom[cx:cx+3]))
                if atom[ct-3] in AA3letter or notProt:
                    self.Write(fmt.format('ATOM  ',iatom,atom[ct-1],atom[ct-3].strip(),    \
                        atom[ct-2].strip(),atom[ct-4].rjust(4),xyz[0],xyz[1],xyz[2],atom[cx+3], \
                        Biso,atom[ct].rjust(12)))
                    if atom[ct-1] == 'OXT':
                        iatom += 1
                        self.Write('{:6s}{:5d}  {:4s}{:3s}'.format('TER   ',iatom,atom[ct-1],atom[ct-3].strip()))
                        nTer += 1
                else:
                    nHet += 1
                    self.Write(fmt.format('HETATM',iatom,atom[ct-1],atom[ct-3].strip(),    \
                        atom[ct-2].strip(),atom[ct-4].rjust(4),xyz[0],xyz[1],xyz[2],atom[cx+3], \
                        Biso,atom[ct].rjust(12)))
                #if atim[cia] == 'a':
                #   put in 'ANISOU' record
                #'ANISOU    1  N   ALA A 340     4392   4159   4615    249   -189     73       N'  
                iatom += 1
            
            vals = [3,0,nHet,0,0,0,6,len(Atoms),nTer,0,nSeq]
            fmt = 'MASTER'+11*'{:5d}'
            self.Write(fmt.format(*vals))
            self.Write('END')
            self.CloseFile()
            print('Phase '+phasenam+' written to PDB file '+self.fullpath)

class ExportPhaseCartXYZ(G2fil.ExportBaseclass):
    '''Used to create a Cartesian XYZ file for a phase

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'Cartesian XYZ',
            extension='.XYZ',
            longFormatName = 'Export phase with Cartesian coordinates as .XYZ file'
            )
        self.exporttype = ['phase']
        self.multiple = True

    def Exporter(self,event=None):
        '''Export as a XYZ file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()
        if self.ExportSelect():    # set export parameters; ask for file name
            return
        filename = self.filename
        self.OpenFile()
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info
            General = phasedict['General']
            i = self.Phases[phasenam]['pId']
            Atoms = phasedict['Atoms']
            if not len(Atoms):
                print('**** ERROR - Phase '+phasenam+' has no atoms! ****')
                continue
            if len(self.phasenam) > 1: # if more than one filename is included, add a phase #
                self.filename = os.path.splitext(filename)[1] + "_" + str(i) + self.extension
            cx,ct,cs,cia = General['AtomPtrs']
            Cell = General['Cell'][1:7]
            A,B = G2lat.cell2AB(Cell)
            fmt = '{:4s}'+3*'{:12.4f}'
            self.Write('{:6d}'.format(len(Atoms)))
            self.Write(' ')
            for atom in Atoms:
                xyz = np.inner(A,np.array(atom[cx:cx+3]))
                self.Write(fmt.format(atom[ct],*xyz))
            self.CloseFile()
            print('Phase '+phasenam+' written to XYZ file '+self.fullpath)
            
class ExportDrawPhaseCartXYZ(G2fil.ExportBaseclass):
    '''Used to create a Cartesian XYZ file for a phase draw atoms

    :param wx.Frame G2frame: reference to main GSAS-II frame
    '''
    def __init__(self,G2frame):
        super(self.__class__,self).__init__( # fancy way to say <parentclass>.__init__
            G2frame=G2frame,
            formatName = 'Draw Atoms Cartesian XYZ',
            extension='.XYZ',
            longFormatName = 'Export draw atoms with Cartesian coordinates as .XYZ file'
            )
        self.exporttype = ['phase']
        self.multiple = True

    def Exporter(self,event=None):
        
        def getRadius(atom):
            atNum = General['AtomTypes'].index(atom[ct])
            if 'H' == atom[ct]:
                if drawingData['showHydrogen']:
                    if 'vdW' in atom[cs]:
                        radius = vdwScale*vdWRadii[atNum]
                    else:
                        radius = ballScale*drawingData['sizeH']
                else:
                    radius = 0.0
            else:
                if 'vdW' in atom[cs]:
                    radius = vdwScale*vdWRadii[atNum]
                else:
                    radius = ballScale*BondRadii[atNum]
            return radius
            
        '''Export as a XYZ file
        '''
        # the export process starts here
        self.InitExport(event)
        # load all of the tree into a set of dicts
        self.loadTree()
        # create a dict with refined values and their uncertainties
        self.loadParmDict()
        if self.ExportSelect():    # set export parameters; ask for file name
            return
        filename = self.filename
        self.OpenFile()
        for phasenam in self.phasenam:
            phasedict = self.Phases[phasenam] # pointer to current phase info
            General = phasedict['General']
            i = self.Phases[phasenam]['pId']
            drawingData = phasedict['Drawing']
            vdwScale = drawingData['vdwScale']
            vdWRadii = General['vdWRadii']
            BondRadii = General['BondRadii']
            bondR = drawingData['bondRadius']
            ballScale = drawingData['ballScale']
            cx,ct,cs,ci = drawingData['atomPtrs']
            Atoms = drawingData['Atoms']
            if not len(Atoms):
                print('**** ERROR - Phase '+phasenam+' has no atoms! ****')
                continue
            if len(self.phasenam) > 1: # if more than one filename is included, add a phase #
                self.filename = os.path.splitext(filename)[1] + "_" + str(i) + self.extension
            Cell = General['Cell'][1:7]
            A,B = G2lat.cell2AB(Cell)
            fmt = '{:4s}'+4*'{:10.3f}'
            line = ' line '+6*'{:10.3f}'
            stick = ' stick '+7*'{:10.3f}'
            face = ' tri '+9*'{:10.3f}'
            self.Write('XYZ file of drawn %s - Cartesian axes suitable for 3D printing/glass etching'%phasenam)
            self.Write('Atoms as balls: element X Y Z radius')
            for atom in Atoms:
                radius = getRadius(atom)
                xyz = np.inner(A,np.array(atom[cx:cx+3]))
                self.Write(fmt.format(atom[ct],*xyz,radius))
            if drawingData['unitCellBox']:
                uBox = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
                uEdges = np.array([
                    [uBox[0],uBox[1]],[uBox[0],uBox[3]],[uBox[0],uBox[4]],[uBox[1],uBox[2]], 
                    [uBox[2],uBox[3]],[uBox[1],uBox[5]],[uBox[2],uBox[6]],[uBox[3],uBox[7]], 
                    [uBox[4],uBox[5]],[uBox[5],uBox[6]],[uBox[6],uBox[7]],[uBox[7],uBox[4]]])
                self.Write('cell edges:')
                xyz = [[0,0,0],[1,1,1]]
                for edge in uEdges:
                    xyz = [np.inner(A,edge[0]),np.inner(A,edge[1])]
                    self.Write(line.format(xyz[0][0],xyz[0][1],xyz[0][2],xyz[1][0],xyz[1][1],xyz[1][2]))
            self.Write('bonds as ball surface to midpoint: x0 y0 z0 x1 y1 z1 r')
            for atom in Atoms:
                xyz = np.inner(A,np.array(atom[cx:cx+3]))
                radius = getRadius(atom)
                Bonds = atom[-2]
                for bond in Bonds:
                    vec = np.inner(A,bond)
                    dist = np.sqrt(np.sum(vec**2))
                    xyz0 = xyz+vec*(radius/dist)
                    bxyz = xyz+vec
                    if 'sticks' in atom[cs]:
                        self.Write(stick.format(xyz0[0],xyz0[1],xyz0[2],bxyz[0],bxyz[1],bxyz[2],bondR))
            self.Write('polygon faces as triangles: x0 y0 z0 x1 y1 z1 x2 y2 z2')
            for atom in Atoms:
                xyz = np.inner(A,np.array(atom[cx:cx+3]))
                Faces = atom[-1]
                for facet in Faces:
                    vx = facet[0]+xyz
                    self.Write(face.format(vx[0,0],vx[0,1],vx[0,2],vx[1,0],vx[1,1],vx[1,2],vx[2,0],vx[2,1],vx[2,2]))
            self.CloseFile()
            print('Phase Draw Atoms '+phasenam+' written to XYZ file '+self.fullpath)
            
                
