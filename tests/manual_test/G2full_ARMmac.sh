WORKSPACE=/tmp
condaHome=/tmp/conda311

gitInstallRepo=git@github.com:AdvancedPhotonSource/GSAS-II-buildtools.git
gitCodeRepo=https://github.com/AdvancedPhotonSource/GSAS-II.git

pyver=3.11
numpyver=1.26

packages="python=$pyver wxpython numpy=$numpyver scipy matplotlib pyopengl conda anaconda-client git gitpython requests pillow h5py imageio scons cython seekpath"

env=testG2
miniforge=https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh

install=True
gitinstall=True

if [ "$install" = "True" ]; then
	rm -rf $condaHome
	if [ ! -e "/tmp/Miniforge-latest.sh" ]; then
		echo Downloading
		curl -L $miniforge -o /tmp/Miniforge-latest.sh
	else
		echo "skipping miniconda download"
	fi
	if [ ! -d "$condaHome" ]; then
		echo creating conda installation
		bash /tmp/Miniforge-latest.sh -b -p $condaHome
	else
		echo "skip miniconda install"
	fi
fi

echo source $condaHome/bin/activate
source $condaHome/bin/activate

if [ "$install" = "True" ]; then
	conda create -y -n $env $packages
fi

set +x
echo source $condaHome/bin/activate $env
source $condaHome/bin/activate $env

if [ "$gitinstall" = "True" ]; then
	set -x
	rm -rf $WORKSPACE/GSAS2-code
	mkdir $WORKSPACE/GSAS2-code
	git clone -b develop $gitCodeRepo $WORKSPACE/GSAS2-code --depth 50
fi

directory="$WORKSPACE/GSAS2-code/GSASII/bindist"

if [ ! -d "$directory" ]; then
	mkdir "$directory"
	echo "Directory $directory created."
else
	echo "Directory $directory already exists."
fi

download_binaries() {
	$condaHome/bin/python - <<END
def download_binaries():
    import urllib.request
    import tarfile
    import os
    
    thetarfile = "https://powder.ornl.gov/files/mac_arm_p3.11_n1.26.tgz"
    
    os.chdir("${WORKSPACE}/GSAS2-code/GSASII/bindist")
    
    ftpstream = urllib.request.urlopen(thetarfile)
    thetarfile = tarfile.open(fileobj=ftpstream, mode="r|gz")
    thetarfile.extractall()

download_binaries()
END
}

download_binaries

python $WORKSPACE/GSAS2-code/GSASII/GSASII.py

exit
