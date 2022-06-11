import re
import copy
import GSASIIobj as G2obj

# compile regular expressions
parenRE = re.compile(r'\)|\(')
bracketRE = re.compile(r'\[|\]')
# get parameter table
G2obj.CompileVarDesc()
fp = open('docs/source/vars.rst','w')
fp.write('''.. 
    This file is created using the makeVarTbl.py. Edit that not this file.
         
.. list-table:: Naming for GSAS-II parameter names, ``p:h:<var>:n``
   :widths: 35 65
   :header-rows: 1

   * - ``<var>``
     - usage
''')

explain = {71:'one or more digits (0, 1,... 9)',
           72:'the digits 0, 1, or 2',
           73:'a digit between 0 and 5'}
regList = {'0-9':71, '0-2':72, '0-5':73}
exmplDig = {71:'10', 72:'0', 73:'0'}
nextChar = 75 # N.B. 74 saved for term at end of line (.*)

for line,r in enumerate(G2obj.reVarDesc): # loop over each entry in table
    comment = ''
    parmName = orig_str = str(r).split("'")[1]
    parmName = parmName.replace(r'\\(','(').replace(r'\\)',')')
    parmName = parmName.replace('[0-9]*','[0-9]')
    if parmName.endswith('$'): # final character in RE
        parmName = parmName.replace('$','')
    termcount = parmName.count('[')
    if parmName.count('(') != 0 and termcount > 0:  # () terms replace \\1, \\2,...
        symTerms = []
        exmplTerms = []
        for i,field in enumerate(parenRE.split(parmName)[1::2]):
            symList = []
            exmplList = []            
            for j,reg in enumerate(bracketRE.split(field)[1::2]):
                if reg not in regList.keys():
                    s = 'one of the characters'
                    for c in reg: 
                        if c == reg[-1]:
                            s = s[:-1] + ' or ' + c
                        else:
                            s += ' ' + c + ','
                    explain[nextChar] = s
                    regList[reg] = nextChar
                    exmplDig[nextChar] = reg[0]
                    nextChar += 1
                if termcount > 1:
                    sym = '\\ :math:`\\scriptstyle '+chr(regList[reg])+'_'+str(j)+'`\\ '
                else:
                    sym = '\\ :math:`\\scriptstyle '+chr(regList[reg])+'`\\ '
                symList.append(sym) 
                exmplList.append(exmplDig[regList[reg]])
                if comment: comment += " and "
                comment += sym + ' is ' + explain[regList[reg]]
            repTmp = copy.copy(bracketRE.split(field))
            repTmp[1::2] = symList
            symTerms.append(''.join(repTmp))
            exmplTmp = copy.copy(bracketRE.split(field))
            exmplTmp[1::2] = exmplList
            exmplTerms.append(''.join(exmplTmp))
            #G2obj.reVarDesc[r].replace()

        repTmp = copy.copy(parenRE.split(parmName))
        repTmp[1::2] = symTerms
        repTmp = ''.join(repTmp)
        exmplTmp = copy.copy(parenRE.split(parmName))
        exmplTmp[1::2] = exmplTerms
        exmplTmp = ''.join(exmplTmp)
        #print(parmName,G2obj.reVarDesc[r])
        out1 = repTmp + ' (example: ``'+exmplTmp+'``)'
        out2 = G2obj.reVarDesc[r].replace('\\'+str(i+1),symTerms[i])
        if comment:
            out2 += "; where " + comment + ','
    elif parmName.endswith('(.*)') or parmName.endswith('.*'):
        sym = '\\ :math:`\\scriptstyle '+chr(74)+'`\\ '
        out1 = repTmp + ' (example: ``'+exmplTmp+'``)'
        comment += ' a number or string ('+sym+') is appended after the semicolon'
        parmName = parmName.replace('(.*)','').replace('.*','')
        out1 = parmName + sym +' (example: ``'+parmName+'11``)'
        out2 = G2obj.reVarDesc[r].replace('\\1',sym)
        if parmName.startswith('Bk'):
            out2 += '; where '+sym+' is the background peak number'
        elif parmName.startswith('Back'):
            out2 += '; where '+sym+' is the background term number'
        elif parmName.startswith('Mustr'):
            out2 += ('; where '+sym+' can be i for isotropic or equatorial'+
                     ' and a is axial (uniaxial broadening),'+
                     ' a number for generalized (Stephens) broadening'+
                     ' or mx for the Gaussian/Lorentzian mixing term (LGmix)')
        elif parmName.startswith('Size'):
            out2 += ('; where '+sym+' can be i for isotropic or equatorial,'+
                     ' and a is axial (uniaxial broadening),'+
                     ' a number between 0 and 5 for ellipsoidal broadening'+
                     ' or mx for the Gaussian/Lorentzian mixing term (LGmix)')
    else:
        out1 = parmName 
        out2 = G2obj.reVarDesc[r]
    if out2[-1] == ',': 
        out2 = out2[:-1] + '.'
    elif out2[-1] != '.':
        out2 += '.'
        
    fp.write('   * - '+out1+'\n')
    fp.write('     - '+out2+'\n')
fp.close()
