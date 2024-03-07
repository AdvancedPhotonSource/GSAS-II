#!/usr/bin/env python
'''
*macRelink: Remove hardcoded library references*
===================================================

'''
from __future__ import division, print_function
import sys, os, glob, subprocess
#os.path, stat, shutil, subprocess, plistlib
def Usage():
    print("\n\tUsage: python "+sys.argv[0]+" <binary-dir>\n")
    sys.exit()
    
def MakeByte2str(arg):
    '''Convert output from subprocess pipes (bytes) to str (unicode) in Python 3.
    In Python 2: Leaves output alone (already str). 
    Leaves stuff of other types alone (including unicode in Py2)
    Works recursively for string-like stuff in nested loops and tuples.

    typical use::

        out = MakeByte2str(out)

    or::

        out,err = MakeByte2str(s.communicate())
    
    '''
    if isinstance(arg,str): return arg
    if isinstance(arg,bytes): return arg.decode()
    if isinstance(arg,list):
        return [MakeByte2str(i) for i in arg]
    if isinstance(arg,tuple):
        return tuple([MakeByte2str(i) for i in arg])
    return arg

if __name__ == '__main__':
    if len(sys.argv) == 2:
        dirloc = os.path.abspath(sys.argv[1])
    else:
        Usage()
        raise Exception

    fileList = glob.glob(os.path.join(dirloc,'*.so'))
    fileList += glob.glob(os.path.join(dirloc,'convcell'))
    fileList += glob.glob(os.path.join(dirloc,'LATTIC'))
    
    libs = {}
    ignorelist = []
    for f in fileList:
        cmd = ['otool','-L',f]
        s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = MakeByte2str(s.communicate())
        for i in out.split('\n')[1:]: 
            if not i: continue
            lib = i.split()[0]
            if os.path.split(lib)[0] == "/usr/lib": # ignore items in system library (usually libSystem.B.dylib)
                if "libSystem" not in lib: print("ignoring ",lib)
                continue
            if "@rpath" in lib and lib not in ignorelist:
                ignorelist.append(lib)
                print("ignoring ",lib)
                continue
            if lib not in libs:
                libs[lib] = []
            libs[lib].append(f)
    for key in libs:
        newkey = os.path.join('@rpath',os.path.split(key)[1])
        print('Fixing',key,'to',newkey)
        for f in libs[key]:
            print('\t',os.path.split(f)[1])
            cmd = ["install_name_tool","-change",key,newkey,f]
            #print(cmd)
            s = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = MakeByte2str(s.communicate())
