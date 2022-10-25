#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
'''
*Module nistlat: NIST*LATTICE cell computations*
------------------------------------------------------

This implements an interface to the NIST*LATTICE code using 
the Spring 1991 program version. NIST*LATTICE, "A Program to Analyze 
Lattice Relationships" was created by Vicky Lynn Karen and Alan D. Mighell, 
National Institute of Standards and Technology, Materials Science and 
Engineering Laboratory, Gaithersburg, Maryland 20899. 

[cite: V. L. Karen and A. D. Mighell, NIST Technical Note 1290 
(1991), https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1290.pdf;
V. L. Karen & A. D. Mighell, U.S. Patent 5,235,523,
https://patents.google.com/patent/US5235523A/en?oq=5235523]

This is still under development; not yet in use.

Easy part done: cell reduction and reversion of that. Work to come, 
"Search for higher symmetry unit cells" & "Relate two unit cells"
'''

import subprocess
import re
import numpy as np

nistlattice = "/Users/toby/boxGSASII/bin/LATTIC"
convcell = "/Users/toby/boxGSASII/bin/convcell"

centerLbl = {'P':'Primitive', 'A':'A-centered', 'B':'B-centered', 'C':'C-centered',
	'F':'Face-centered', 'I':' Body-centered', 'R':'Rhombohedral', 'H':'Rhombohedral'}

def showCell(cell,center,setting):
    'format cell input or output'
    s = "{:.4f} {:.4f} {:.4f} {:.3f} {:.3f} {:.3f}".format(*cell)
    s += '  ' + centerLbl.get(center,'?')
    if setting.strip(): s += '-' + setting
    return s

def uniqCells(cellList):
    'remove duplicated cells from an output list'
    uList = []
    shown = []
    for i in cellList:
        s = "{:} {:.3f} {:.3f} {:.3f} {:.2f} {:.2f} {:.2f}".format(i[4],*i[1])
        if s in shown: continue
        shown.append(s)
        uList.append(i)
    return uList
        
def ReduceCell(center,cellin,mode=0,deltaV=0):
    '''Compute reduced cell(s) with NIST*LATTICE
    
    :param str center: cell centering code, one of P/A/B/C/F/I/R
      Note that 'R' is used for rhombohedral lattices in either 
      hexagonal or rhombohedral (primitive) cells
    :param list cellin: six lattice constants as float values
    :param int mode: 
        0: reduction, 
        1: generate supercells, 
        2: generate subcells
        3: generate sub- and supercells
    :param int deltaV: volume ratios for sub/supercells if mode != 0 as 
      ratio of original cell to smallest subcell or largest supercell. Ignored 
      if mode=0. 
    :returns: a dict with item 'input' with input cell as (cell,center,setting)
      and 'output' which is a list of reduced cells of form
      (d,cell,vol,mat,center,setting). In these, 
        cell: is the six cell dimensions;
        center: is as above (always 'P' on output); 
        setting: is ' ' except for rhombohedral symmetry where it may be R or H for the cell type;
        d: is the volume ratio for new cell over input cell;
        vol: is volume of output cell
        mat: is the matrix that gives the output cell when the input cell is multiplied by mat.     
    '''

    # setting is non-blank for rhombohedral lattices (center=R) only.
    #  setting = R for rhob axes 
    #  setting = H for hex. axes
    setting = " " 
    (a,b,c,alpha,beta,gamma) = cellin
    celldict = {}
    if center == "R":
        setting = "R"
    elif alpha == 90 and beta == 90 and gamma == 120:
        setting = "H"
    # prepare input and start program
    cellline = '{:10.4f}{:10.4f}{:10.4f}{:10.3f}{:10.3f}{:10.3f}'.format(*cellin)
    inp = "RSS      1\n"
    inp += "{:1d}  {:1d}{:1d}   {}{}{}\n".format(mode,2,deltaV,center,setting,cellline)
    p = subprocess.Popen([nistlattice],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p.stdin.write(bytearray(inp,'utf8'))
    p.stdin.close()
    # read output and parse
    err = p.stderr.read()
    celldict['input'] = (cellin,center,setting)
    celldict['output'] = []
    d = 1
    line = '?'
    linenum = 0
    try:
        for b in p.stdout.readlines():
            linenum += 1
            line = b.decode()
            pat = r"T 2= (.*)/ (.*)/ (.*)"  # transform matrix
            s = re.split(pat,line)
            if len(s) > 1: 
                mat = [[float(j) for j in re.split(r" *([\d\.-]*) *",i,maxsplit=3)[1::2]]
                           for i in s[1:-1]]
            
            pat = r"Delta =(.*)"   # Volume ratio (if mode >0)
            s = re.split(pat,line)
            if len(s) > 1: 
                d = float(eval(s[1]))

            pat = r"CELL  2=(.*)V2=(.*)"  # cell lattice and volume
            s = re.split(pat,line)
            if len(s) > 1:
                lat = [float(i) for i in re.split(r" *([\d\.-]*) *",s[1],maxsplit=6)[1::2]]
                vol = float(re.split(r" *([\d\.-]*) *",s[2],maxsplit=1)[1])
                celldict['output'].append((d,lat,vol,mat,'P',' ')) # note reduced cells are all primitive            
    except:
        print('Parse error at line ',linenum)
        print(line)
        return celldict
    finally:
        p.terminate()
    if len(celldict['output'] or len(err) > 0) == 0:
        print('Error:' ,err.decode())
    return celldict

def ConvCell(redcell):
    '''Converts a reduced cell to a conventional cell
    :param list redcell: unit cell parameters as 3 cell lengths 
      and 3 angles (in degrees)
    :returns: tuple (cell,center,setting,mat) where 
        cell: has the six cell dimensions for the conventional cell;
        center: is P/A/B/C/F/I/R;
        setting: is ' ' except for rhombohedral symmetry (center=R), where 
          it will always be H (for hexagonal cell choice);
        mat: is the matrix that gives the conventional cell when the reduced 
          cell is multiplied by mat.
    '''
    inp = "{:10.5f}{:10.5f}{:10.5f}{:10.4f}{:10.4f}{:10.4f}".format(*redcell)
    p = subprocess.Popen([convcell],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p.stdin.write(bytearray(inp,'utf8'))
    p.stdin.close()
    # read output and parse
    err = p.stderr.read()
    line = '?'
    linenum = 0
    cell = []
    center = ' '
    setting = ' '
    try:
        for b in p.stdout.readlines():
            linenum += 1
            line = b.decode()
            if linenum == 1:
                cell = [float(i) for i in line.split()[:6]]
                center = line.split()[-1]
                if center == 'H':
                    center = 'R'
                    setting = 'H'
            if linenum == 2:
                mat = np.array([float(i) for i in line.split()]).reshape(3,3)
    except:
        print('Parse error at line ',linenum)
        print(line)
        #return cell
    finally:
        p.terminate()
    if len(err) >  0:
        print('Error:' ,err.decode())
    return (cell,center,setting,mat)


if __name__ == '__main__':  # test code
    import GSASIIlattice as G2lat
    cellin = [5.,5.,5.,85.,85.,85.,]
    cellList = ConvCell(cellin)
    print('Input reduced cell',showCell(cellin,'P',' '),'\nout',showCell(*cellList[:3]))
    mat = cellList[-1]
    print('test with\n',mat)
    print(G2lat.TransformCell(cellin,mat))
    
    print('\ntest from [5,6,7,90,105,90] C-Centered') # ==>  [5,6,7,90,105,90] C-Centered
    cellin = [3.9051,3.9051,7,99.537,99.537,100.389]
    cellList = ConvCell(cellin)
    print('Input reduced cell',showCell(cellin,'P',' '),'\nout',showCell(*cellList[:3]))
    mat = cellList[-1]
    print('test with\n',mat)
    print(G2lat.TransformCell(cellin,mat))
    
    print('\nHexagonal test (no change)')
    cellin = [5.,5.,7.,90.,90.,120.,]
    cellList = ConvCell(cellin)
    print('Input reduced cell',showCell(cellin,'P',' '),'\nout',showCell(*cellList[:3]))
    mat = cellList[-1]
    print('test with\n',mat)
    print(G2lat.TransformCell(cellin,mat))
    
            
    # cell = ReduceCell('F',[3,3,5,90,90,90],1,2)
    # print('supercell',showCell(*cell['input']))
    # for i in cell['output']: print(i)
    # print('subcell',showCell(*cell['input']))
    # cell = ReduceCell('F',[3,3,5,90,90,90],2,2)
    # for i in cell['output']: print(i)
    # print('both',showCell(*cell['input']))
    cell = ReduceCell('F',[3,3,5,90,90,90],3,2)
    print('\nunique from both mode',showCell(*cell['input']))
    for i in uniqCells(cell['output']): print(i)
    
    cell = ReduceCell('I',[3,3,5,90,90,90])
    print('\ndefault',showCell(*cell['input']))
    for i in cell['output']: print(i)
    print('invert back',ConvCell(cell['output'][0][1]))
    
    
    print('\nF-tetragonal is not standard setting')
    cell = ReduceCell('F',[3,3,5,90,90,90])
    print('default',showCell(*cell['input']))
    for i in cell['output']: print(i)
    print('invert back',ConvCell(cell['output'][0][1]))
