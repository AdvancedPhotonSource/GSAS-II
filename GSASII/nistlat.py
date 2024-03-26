#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date: 2023-11-14 10:44:24 -0600 (Tue, 14 Nov 2023) $
# $Author: toby $
# $Revision: 5701 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/nistlat.py $
# $Id: nistlat.py 5701 2023-11-14 16:44:24Z toby $
########### SVN repository information ###################
'''
This implements an interface to the NIST*LATTICE code using 
the Spring 1991 program version. NIST*LATTICE, "A Program to Analyze 
Lattice Relationships" was created by Vicky Lynn Karen and Alan D. Mighell
(National Institute of Standards and Technology, Materials Science and 
Engineering Laboratory, Gaithersburg, Maryland 20899.)
Minor code modifications made to provide more significant digits for
cell reduction matrix terms. 

Please cite V. L. Karen and A. D. Mighell, NIST Technical Note 1290 (1991),
https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1290.pdf;
and V. L. Karen & A. D. Mighell, U.S. Patent 5,235,523,
https://patents.google.com/patent/US5235523A/en?oq=5235523 if this module 
is used. 

'''

import os
import os.path
import subprocess
import re
import numpy as np
import GSASIIpath
GSASIIpath.SetBinaryPath()
import GSASIIlattice as G2lat

nistlattice = os.path.join(GSASIIpath.binaryPath,"LATTIC")
convcell = os.path.join(GSASIIpath.binaryPath,"convcell")
is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)

debug = False
#debug = True

centerLbl = {'P':'Primitive', 'A':'A-centered', 'B':'B-centered', 'C':'C-centered',
	'F':'Face-centered', 'I':' Body-centered', 'R':'Rhombohedral', 'H':'Rhombohedral'}
   
def showCell(cell,center='P',setting=' ',*ignored):
    '''show unit cell input or output nicely formatted.

    :param list cell: six lattice constants as float values; a 7th volume value
      is ignored if present.
    :param str center: cell centering code; one of P/A/B/C/F/I/R
      Note that 'R' is used for rhombohedral lattices in either 
      rhombohedral (primitive) or hexagonal cells.
    :param str setting: is ' ' except for rhombohedral symmetry where 
      it will be R or H for the cell type.
    :returns: a formatted string
    '''
    s = "{:.4f} {:.4f} {:.4f} {:.3f} {:.3f} {:.3f}".format(*cell[:6])
    s += '  ' + centerLbl.get(center,'?')
    if setting.strip(): s += '-' + setting
    return s

def printCell(label,*args,**kwargs):
    print(label,showCell(*args,**kwargs))

def uniqCells(cellList):
    '''remove duplicated cells from a cell output list from :func:`ReduceCell`

    :param list cellList: A list of reduced cells where each entry represents a
      reduced cell as (_,cell,_,_,center,...) where cell has six lattice 
      constants and center is the cell centering code (P/A/B/C/F/I/R).
    :returns: a list as above, but where each unique cell is listed only once
    '''
    uList = []
    shown = []
    for i in cellList:
        s = "{:} {:.3f} {:.3f} {:.3f} {:.2f} {:.2f} {:.2f}".format(i[4],*i[1])
        if s in shown: continue
        shown.append(s)
        uList.append(i)
    return uList

def _emulateLP(line,fp):
    '''Emulate an antique 132 column line printer, where the first column
    that is printed is used for printer control. '1' starts a new page 
    and '0' double-spaces. Not implemented is '+' which overprints.
    If the file pointer is not None, the line is copied to the 
    file, removing the first column but adding a bit of extra formatting 
    to separate pages or add a blank line where needed. 

    :param str line: string to be printed
    :param file fp: file pointer object
    '''
    if not fp: return
    if line[0] == '1':
        fp.write('\n'+130*'='+'\n'+line[1:])            
    elif line[0] == '0':
        fp.write('\n'+line[1:])            
    elif len(line) == 1:
        fp.write('\n')            
    else:
        fp.write(line[1:])            
        
def ReduceCell(center, cellin, mode=0, deltaV=0, output=None):
    '''Compute reduced cell(s) with NIST*LATTICE
    
    :param str center: cell centering code; one of P/A/B/C/F/I/R
      Note that 'R' is used for rhombohedral lattices in either 
      hexagonal or rhombohedral (primitive) cells
    :param list cellin: six lattice constants as float values
    :param int mode: 

        * 0: reduction, 
        * 1: generate supercells, 
        * 2: generate subcells
        * 3: generate sub- and supercells

    :param int deltaV: volume ratios for sub/supercells if mode != 0 as 
      ratio of original cell to smallest subcell or largest supercell 
      to original cell. Ignored if mode=0. Otherwise should be 2, 3, 4 or 5
    :param str output: name of file to write the NIST*LATTICE output.
      Default is None, which does not produce a file.    
    :returns: a dict with two items, 'input' and 'output'. The value for
      'input' is the input cell as (cell,center,setting). The value for
      'output' is a list of reduced cells of form
      (d,cell,vol,mat,center,setting). In these: 

        * cell: a list with the six cell dimensions;
        * center: is as above (always 'P' on output); 
        * setting: is ' ' except for rhombohedral symmetry where it may be R or H for the cell type;
        * d: is the volume ratio for new cell over input cell;
        * vol: is volume of output cell
        * mat: is the matrix that gives the output cell when the input cell is multiplied by mat.
    '''

    # setting is non-blank for rhombohedral lattices (center=R) only.
    #  setting = R for rhob axes 
    #  setting = H for hex. axes
    setting = " " 
    (a,b,c,alpha,beta,gamma) = cellin
    celldict = {}
    if center == "R" and alpha == 90 and beta == 90 and gamma == 120:
        setting = "H"
    elif center == "R" :
        setting = "R"
    # prepare input and start program
    cellline = '{:10.4f}{:10.4f}{:10.4f}{:10.3f}{:10.3f}{:10.3f}'.format(*cellin)
    inp = "RSS      1\n"
    inp += "{:1d}  {:1d}{:1d}   {}{}{}\n".format(mode,2,deltaV,center,setting,cellline)
    if os.path.exists('NIST10'): # cleanup
        print("Removing old NIST10 file")
        os.remove('NIST10')
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
    fp = None
    if output: fp = open(output,'w')
    try:
        for b in p.stdout.readlines():
            linenum += 1
            line = b.decode()
            _emulateLP(line,fp)
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
        print('ReduceCell parse error at line ',linenum)
        print(line)
        return celldict
    finally:
        p.terminate()
    if len(celldict['output']) == 0 or len(err) > 0:
        print('Error:' ,err.decode())
    return celldict

def ConvCell(redcell):
    '''Converts a reduced cell to a conventional cell

    :param list redcell: unit cell parameters as 3 cell lengths 
      and 3 angles (in degrees)
    :returns: tuple (cell,center,setting,mat) where:
 
        * cell: has the six cell dimensions for the conventional cell;
        * center: is P/A/B/C/F/I/R;
        * setting: is ' ' except for rhombohedral symmetry (center=R), where 
          it will always be H (for hexagonal cell choice);
        * mat: is the matrix that gives the conventional cell when the reduced 
          cell is multiplied by mat.
    '''
    inp = "{:10.5f}{:10.5f}{:10.5f}{:10.4f}{:10.4f}{:10.4f}".format(*redcell)
    if os.path.exists('NIST10'): # cleanup
        print("Removing old NIST10 file")
        os.remove('NIST10')
    p = subprocess.Popen([convcell],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p.stdin.write(bytearray(inp,'utf8'))
    p.stdin.close()
    # read output and parse
    err = p.stderr.read()
    if debug and err:
        print('ConvCell err=',err)
    line = '?'
    linenum = 0
    cell = []
    center = ' '
    setting = ' '
    try:
        for b in p.stdout.readlines():
            line = b.decode()
            if '**WARNING**' in line: 
                print('Note: Warning generated in conversion of reduced\n  cell',
                      redcell,'\n  (Probably OK to ignore)')
                continue
            if not line.strip(): continue
            linenum += 1
            if linenum == 1:
                cell = [float(i) for i in line.split()[:6]]
                center = line.split()[-1]
                if center == 'H':
                    center = 'R'
                    setting = 'H'
            if linenum == 2:
                mat = np.array([float(i) for i in line.split()]).reshape(3,3)
    except:
        print('ConvCell parse error at line ',linenum)
        print(line)
        if debug:
            print("\nRemaining lines:")
            for b1 in p.stdout.readlines():
                print(b1.decode())
        return None
        #return cell
    finally:
        p.terminate()
    if len(err) >  0:
        print('Error:' ,err.decode())
    return (cell,center,setting,mat)


def CompareCell(cell1, center1, cell2, center2, tolerance=3*[0.2]+3*[1],
                    mode='I', vrange=8, output=None):
    '''Search for matrices that relate two unit cells 
    
    :param list cell1: six lattice constants as float values for 1st cell
    :param str center1: cell centering code for 1st cell; one of P/A/B/C/F/I/R
      Note that 'R' is used for rhombohedral lattices in either 
      hexagonal or rhombohedral (primitive) cells
    :param list cell2: six lattice constants as float values for 2nd cell
    :param str center2: cell centering code for 2nd cell (see center1)
    :param list tolerance: comparison tolerances for a, b, c, alpha, beta 
      & gamma (defaults to [0.2,0.2,0.2,1.,1.,1.]
    :param str mode: search mode, which should be either 'I' or 'F'
      'I' provides searching with integral matrices or
      'F' provides searching with integral and fractional matrices
    :param int vrange: maximum matrix term range. 
       Must be 1 <= vrange <= 10 for mode='F' or 
       Must be 1 <= vrange <= 40 for mode='I' 
    :param str output: name of file to write the NIST*LATTICE output.
      Default is None, which does not produce a file.    

    :returns: A list of matrices that match cell1 to cell2 where 
      each entry contains (det, im, m, tol, one2two, two2one) where:

        * det is the determinant, giving the volume ratio between cells
        * im relates the reduced cell for cell1 to the reduced cell for cell2
        * m relates the reduced cell for cell2 to the reduced cell for cell1
        * tol shows the quality of agreement, as six differences between the 
             two reduced cells
        * one2two: a numpy matrix that transforms cell1 to cell2
        * two2one: a numpy matrix that transforms cell2 to cell1
    '''
    # reduce input cells. Save cell, volume and matrix
    rcVmat = [[],[]]
    rcVmat[0] = ReduceCell(center1,cell1)['output'][0][1:4]
    rcVmat[1] = ReduceCell(center2,cell2)['output'][0][1:4]
    reverseCells = False
    if rcVmat[0][1] < rcVmat[1][1]: # largest volume goes first
        reverseCells = True
        rcVmat = list(reversed(rcVmat))
    # prepare input and start program
    inp = "REL      1\n"
    inp += "{:.1s}  {:2d}     {:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}\n".format(
            mode,vrange,*tolerance)
    for i in range(2):
        inp += "          {:10.5f}{:10.5f}{:10.5f}{:10.4f}{:10.4f}{:10.4f}\n".format(
            *rcVmat[i][0])
    if os.path.exists('NIST10'): # cleanup
        print("Removing old NIST10 file")
        os.remove('NIST10')
    p = subprocess.Popen([nistlattice],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p.stdin.write(bytearray(inp,'utf8'))
    p.stdin.close()
    err = p.stderr.read()
    line = '?'
    fp = None
    if output: fp = open(output,'w')
    # read output and parse
    lines = [b.decode() for b in p.stdout.readlines()]
    p.terminate()
    if fp:
        for line in lines: _emulateLP(line,fp)
        fp.close()
    lnum = 0
    while lnum < len(lines):
        if ' H Matrix ' in lines[lnum] and ' Determinant' in lines[lnum]:
            lnum += 2
            break
        lnum += 1
    else:
        print('no header found')
        return []

    xforms = []
    while lnum < len(lines):
        try:
            sl = lines[lnum].split()
            if len(sl) == 10:
                sl = [float(i) for i in sl]
                m = 3*[3*[0]]  # forward matrix
                im = 3*[3*[0]] # inverse matrix
                tol = 6*[0]
                m[0] = sl[0:3]
                im[0] = sl[3:6]
                tol[0:3] = sl[6:9]
                det = sl[9]
            
                lnum += 1
                sl = [float(i) for i in lines[lnum].split()]
                m[1] = sl[0:3]
                im[1] = sl[3:6]
                tol[3:] = sl[6:9]
            
                lnum += 1
                sl = [float(i) for i in lines[lnum].split()]
                m[2] = sl[0:3]
                im[2] = sl[3:6]
            
                xmatI = np.dot(np.dot(np.linalg.inv(rcVmat[0][2]),m),rcVmat[1][2])
                xmat  = np.dot(np.dot(np.linalg.inv(rcVmat[1][2]),im),rcVmat[0][2])
                if reverseCells:
                    one2two, two2one = xmatI, xmat
                    xforms.append((det, m, im, tol, one2two, two2one))
                else:
                    one2two, two2one = xmat, xmatI
                    xforms.append((det, im, m, tol, one2two, two2one))
        except Exception as msg:
            print('error with line',lnum,msg)
            print(lines[lnum])

        lnum += 1
    if len(err) > 0:
        print('Execution error:' ,err.decode())
    return xforms

def CellSymSearch(cellin, center, tolerance=3*[0.2]+3*[1], mode=0,
                      deltaV=2, output=None):
    '''Search for a higher symmetry lattice related to an input unit cell,
    and optionally to the supercells and/or subcells with a specified 
    volume ratio to the input cell. 

    :param list cellin: six lattice constants as float values
    :param str center: cell centering code; one of P/A/B/C/F/I/R
      Note that 'R' is used for rhombohedral lattices in either 
      hexagonal or rhombohedral (primitive) cells
    :param list tolerance: comparison tolerances for a, b, c, alpha, beta 
      & gamma (defaults to [0.2,0.2,0.2,1.,1.,1.]
    :param int mode: 

        * 0: use only input cell,
        * 1: generate supercells, 
        * 2: generate subcells
        * 3: generate sub- and supercells

    :param int deltaV: volume ratios for sub/supercells if mode != 0 as 
      ratio of original cell to smallest subcell or largest supercell 
      to original cell. Ignored if mode=0. Otherwise should be 2, 3, 4 or 5
    :param str output: name of file to write the NIST*LATTICE output. Default
      is None, which does not produce a file.

    :returns: a list of processed cells (only one entry in list when mode=0) 
      where for each cell the the following items are included: 

        * conventional input cell; 
        * reduced input cell; 
        * symmetry-generated conventional cell; 
        * symmetry-generated reduced cell; 
        * matrix to convert sym-generated output cell to input conventional cell
    '''
    # setting is non-blank for rhombohedral lattices (center=R) only.
    #  setting = R for rhob axes 
    #  setting = H for hex. axes
    setting = " " 
    (a,b,c,alpha,beta,gamma) = cellin
    celldict = {}
    if center == "R" and alpha == 90 and beta == 90 and gamma == 120:
        setting = "H"
    elif center == "R" :
        setting = "R"
    # prepare input and start program
    cellline = '{:10.4f}{:10.4f}{:10.4f}{:10.3f}{:10.3f}{:10.3f}'.format(*cellin)
    inp = "SYM      1\n"
    inp += "R{:1d} {:2d}     {:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}\n".format(
            mode,deltaV,*tolerance)
    inp += "{:1d}  {:1d}{:1d}   {}{}{}\n".format(mode,2,deltaV,center,setting,cellline)
    if os.path.exists('NIST10'): # cleanup
        print("Removing old NIST10 file")
        os.remove('NIST10')
    p = subprocess.Popen([nistlattice],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p.stdin.write(bytearray(inp,'utf8'))
    p.stdin.close()
    # read output and parse
    err = p.stderr.read()
    
    d = 1
    lines = [b.decode() for b in p.stdout.readlines()]
    p.terminate()
    fp = None
    if output: fp = open(output,'w')
    if fp:
        for line in lines: _emulateLP(line,fp)
        fp.close()
    lnum = 0
    d = 1
    startCellList = []  # reduced cell(s) as (determinant, xform-matrix, reduced-cell, volume) save for internal use
    symCellList = []
    try:
        # get list of starting cells
        while lnum < len(lines):
            line = lines[lnum]
            lnum += 1
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
                startCellList.append((d,mat,lat,vol))
            if "  **  Symmetry" in line:
                break
        #======================================================================
        # for each input-generated cell, get sets of generated H-matrices 
        for icell,cellstuff in enumerate(startCellList):
            while lnum < len(lines):
                line = lines[lnum]
                lnum += 1
                pat = r"CELL  1=(.*)V1=(.*)"  # cell lattice and volume
                s = re.split(pat,line)
                if len(s) > 1:
                    lat = [float(i) for i in re.split(r" *([\d\.-]*) *",s[1],maxsplit=6)[1::2]]
                    vol = float(re.split(r" *([\d\.-]*) *",s[2],maxsplit=1)[1])
                    if not np.allclose(cellstuff[2],lat):
                        print('Warning mismatch between list of reduced cells and'
                                  ' processed output at line',
                                  lnum)
                        print('cell info=',cellstuff)
                        print('lat,vol=',lat,vol)
                    break
            while lnum < len(lines):
                line = lines[lnum]
                lnum += 1
                if ' H Matrix ' in lines[lnum] and ' Determinant' in lines[lnum]:
                    lnum += 2
                    break
            # for each H-matrix, compute a new reduced cell
            xformSum = np.zeros(7)
            xformCount = 0
            while lnum < len(lines):
                line = lines[lnum]
                sl = lines[lnum].split()
                if len(sl) == 10:
                    sl = [float(i) for i in sl]
                    m = 3*[3*[0]]  # forward matrix
                    im = 3*[3*[0]] # inverse matrix
                    tol = 6*[0]
                    m[0] = sl[0:3]
                    im[0] = sl[3:6]
                    tol[0:3] = sl[6:9]
                    #det = sl[9]
            
                    lnum += 1
                    sl = [float(i) for i in lines[lnum].split()]
                    m[1] = sl[0:3]
                    im[1] = sl[3:6]
                    tol[3:] = sl[6:9]
            
                    lnum += 1
                    sl = [float(i) for i in lines[lnum].split()]
                    m[2] = sl[0:3]
                    im[2] = sl[3:6]
                    xformSum += G2lat.TransformCell(lat,m)
                    xformCount += 1
                lnum += 1
                    
                # got to end of output for current cell, process this input    
                if 50*'-' in line:
                    if xformCount:
                        inRedCell = (startCellList[icell][2], 'P', ' ')
                        if icell ==0:
                            inCnvCell = (cellin, center, setting)
                            red2convInp = np.linalg.inv(startCellList[0][1])
                        else:
                            inCnvCell = ConvCell(startCellList[icell][2])
                            red2convInp = inCnvCell[3]
                        redCell = ([j for j in (xformSum/xformCount)[:6]], 'P', ' ')
                        cnvCell = ConvCell(redCell[0][:6])
                        if cnvCell is None:
                            print('ConvCell failed on ',redCell[0][:6])
                            continue
                        # steps to walk back conventional cell back to original
                        #c1 = G2lat.TransformCell(cnvCell[0],np.linalg.inv(cnvCell[3]))
                        #c2 = G2lat.TransformCell(c1[:6],im) # use last, as they are all equiv.
                        #c3 = G2lat.TransformCell(c2[:6],red2convInp)
                        cnv2origMat = np.dot(np.dot(red2convInp,im),np.linalg.inv(cnvCell[3]))
                        #c3 = G2lat.TransformCell(cnvCell[0],cnv2origMat)
                        symCellList.append((
                            inCnvCell[:3],   # input cell (conventional)
                            inRedCell,       # input cell (reduced)
                            cnvCell[:3],   # symmetry-generated conventional cell
                            redCell,       # symmetry-generated reduced cell
                            cnv2origMat,    # matrix to convert sym-generated output cell to input conventional cell
                            ))                            
                    break   # go on to next input-generated cell

        
    except Exception as msg:
        print('CellSymSearch parse error at line ',lnum,'\nNote error:',msg)
        print(line)
        return celldict
    finally:
        p.terminate()
    if len(symCellList) == 0 or len(err) > 0:
        print('Error:' ,err.decode())
    return symCellList
    
if __name__ == '__main__':  # test code

    cell = [14.259, 22.539, 8.741, 90., 114.1, 90.]
    center = 'C'
    print('CellSymSearch output')
    for i in CellSymSearch(cell, center, output='/tmp/cellsym.txt'):
        print('\ninp-conv, inp-red, out-conv, out-red, out(conv)->inp(conv)')
        for j in i: print(j)
            
    print('\n\nCellSymSearch output')
    for i in CellSymSearch(cell, center, output='/tmp/cellsym.txt', mode=3):
        print('\ninp-conv, inp-red, out-conv, out-red, out(conv)->inp(conv)')
        for j in i: print(j)

    import sys; sys.exit()

    cell1 = (5.03461,5.03461,13.74753,90,90,120)
    center1 = 'R'
    cell2 = (7.40242,5.03461,5.42665,90,84.14,90)
    center2 = 'P'
    #print(ReduceCell(center1,cell1,output='/tmp/reduce.txt'))
    tolerance = 3*[0.1]+3*[.5]
    print('\ncomparing ',showCell(cell1,center1,' '),showCell(cell2,center2,' '))
    out = CompareCell(cell1, center1, cell2, center2, tolerance, output='/tmp/CompCell.txt')
    for i in out:
        print(70*'=')
        print(i[0],i[3],'\n',i[5])
    print('testing with 1st entry')
    print('cell1 ',showCell(cell1,center1,' '))
    print('cell2-->cell1',G2lat.TransformCell(cell2,out[0][5]))
    print('cell2 ',showCell(cell2,center2,' '))
    print('cell1-->cell2',G2lat.TransformCell(cell1,out[0][4]))

    out = CompareCell(cell2, center2, cell1, center1, tolerance)
    print('reversed input, testing with 1st entry')
    print('cell1 ',showCell(cell2,center2,' '))
    print('cell2-->cell1',G2lat.TransformCell(cell1,out[0][5]))
    print('cell2 ',showCell(cell1,center1,' '))
    print('cell1-->cell2',G2lat.TransformCell(cell2,out[0][4]))

    
    center2 = 'I'
    print('\n\ncomparing ',showCell(cell1,center1,' '),showCell(cell2,center2,' '))
    out = CompareCell(cell1, center1, cell2, center2, tolerance)
    for i in out:
        print(70*'=')
        print(i[0],i[3],'\n',i[5])
    print('testing with 1st entry')
    print('cell1 ',showCell(cell1,center1,' '))
    print('cell2-->cell1',G2lat.TransformCell(cell2,out[0][5]))
    print('cell2 ',showCell(cell2,center2,' '))
    print('cell1-->cell2',G2lat.TransformCell(cell1,out[0][4]))

    import sys; sys.exit()

    
    cellin = [5.,5.,5.,85.,85.,85.,]
    cellList = ConvCell(cellin)
    if cellList is None: 
        print('ConvCell failed',cellin)
        sys.exit()
    print('Input reduced cell',showCell(cellin,'P',' '),'\nout',showCell(*cellList))
    mat = cellList[-1]
    print('test with\n',mat)
    print(G2lat.TransformCell(cellin,mat))
    
    print('\ntest from [5,6,7,90,105,90] C-Centered') # ==>  [5,6,7,90,105,90] C-Centered
    cellin = [3.9051,3.9051,7,99.537,99.537,100.389]
    cellList = ConvCell(cellin)
    if cellList is None: 
        print('ConvCell failed',cellin)
        sys.exit()
    print('Input reduced cell',showCell(cellin,'P',' '),'\nout',showCell(*cellList))
    mat = cellList[-1]
    print('test with\n',mat)
    print(G2lat.TransformCell(cellin,mat))
    
    print('\nHexagonal test (no change)')
    cellin = [5.,5.,7.,90.,90.,120.,]
    cellList = ConvCell(cellin)
    if cellList is None: 
        print('ConvCell failed',cellin)
        sys.exit()
    print('Input reduced cell',showCell(cellin,'P',' '),'\nout',showCell(*cellList))
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
