# this routine is intended to be run from inside meson. It creates a
# copy of pyspg.for where the latest git version information is inserted
# into the FORTVER subroutine. Usage:
#    python patch_pyspg.py <input pyspg.for> <output file>
# created by copilot (BHT 11/26)
import sys
import os
import re
import shutil
import subprocess
source = os.environ.get('MESON_SOURCE_ROOT',
                    os.path.dirname(os.path.dirname(__file__)))

version = ''
number = '0'
githash = ''
if shutil.which('git'):
    out = subprocess.run(['git','tag','-l','--sort=-authordate','v*'],
                             stdout=subprocess.PIPE,cwd=source)
    version = out.stdout.decode('latin-1').split('\n')[0].strip()
    out = subprocess.run(['git','tag','-l','--sort=-authordate','[0-9]*'],
                             stdout=subprocess.PIPE,cwd=source)
    number = out.stdout.decode('latin-1').split('\n')[0].strip()
    out = subprocess.run(['git','log','HEAD','-1'],
                             stdout=subprocess.PIPE,cwd=source)
    log = out.stdout.decode('latin-1')
    if 'commit' in log:
        githash = log.split('commit')[1].split('\n')[0].strip()
else:
    try:
        sys.path.insert(0,source)
        from GSASII import git_verinfo
        githash = git_verinfo.git_version
        version = git_verinfo.git_versiontag
        for item in git_verinfo.git_tags+git_verinfo.git_prevtags:
            if item.isnumeric():
                number = item
                break
    except:
        pass

if version:    
    #verstr = f"{version}:{githash[:7]}"[:15] # VERSTR is declared as CHARACTER*(15) in pyspg.for
    verstr = f"{version}"
else:
    verstr = "*Undef*"
if not number.isnumeric():
    number = '0'

infile, outfile = sys.argv[1], sys.argv[2]
with open(infile) as fp:
    contents = fp.read()
contents = re.sub(r'(VERSTR\s*=\s*)["\'][^"\']*["\']',
                  rf'\g<1>"{verstr}"', contents)
contents = re.sub(r'VERNUM\s*=\s*\d+', f'VERNUM = {number}', contents)
# don't touch the output file if nothing changed, so that dependent
# build steps are not needlessly rerun
if os.path.exists(outfile):
    with open(outfile) as fp:
        if fp.read() == contents:
            print(f'{outfile} unchanged (version {verstr}, {number})')
            sys.exit()
with open(outfile,'w') as fp:
    fp.write(contents)
print(f'version {verstr} ({number}) written to {outfile}')
