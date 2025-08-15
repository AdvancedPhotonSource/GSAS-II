# this routine is intended to be run from inside meson and creates a
# file named GSASIIversion.txt with that information
import datetime
import os
import shutil
import subprocess
import sys

source = os.environ.get('MESON_SOURCE_ROOT',
                    os.path.dirname(os.path.dirname(__file__)))
build = os.environ.get('MESON_BUILD_ROOT',
                    os.path.dirname(os.path.dirname(__file__)))

if shutil.which('git'):
    out = subprocess.run(['git','tag','-l','--sort=-authordate','v*'],
                             stdout=subprocess.PIPE,cwd=source, check=False)
    version = out.stdout.decode('latin-1').split('\n')[0].strip()
    out = subprocess.run(['git','tag','-l','--sort=-authordate','[0-9]*'],
                             stdout=subprocess.PIPE,cwd=source, check=False)
    number = out.stdout.decode('latin-1').split('\n')[0].strip()
    out = subprocess.run(['git','log','HEAD','-1'],
                             stdout=subprocess.PIPE,cwd=source, check=False)
    githash = out.stdout.decode('latin-1').split('commit')[1].split('\n')[0].strip()
    how = 'git'
else:
    how = 'git_verinfo'
    try:
        sys.path.insert(0,source)
        from GSASII import git_verinfo
        number = '?'
        githash = git_verinfo.git_version
        version = git_verinfo.git_versiontag
        for item in git_verinfo.git_tags+git_verinfo.git_prevtags:
            if item.isnumeric():
                number = item
                break
    except:
        sys.exit()
outfile = os.path.join(build,'sources','GSASIIversion.txt')
with open(outfile,'w') as fp:
    fp.write(f"{version}\n{number}\ncommit {githash}\n")
    fp.write(f"#created by {__file__} on {datetime.datetime.now()} from {how}\n")
    fp.close()
print(outfile)
