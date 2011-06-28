# this file is present for historical reasons. The build is now done with scons, which does not need to 
# be configured as source code files are added or removed
#
scons:
	scons -Q
	@echo "************************************"
	@echo " Note: build is now done with scons"
	@echo "************************************"

# build the compiled fortran codes needed by GSAS-II with f2py
# This is a work in progress. The idea is that files built here are placed in 
# the bin directory, which is empty in the svn distribution. If these files will 
# be distributed, they should be moved to a directory named bin<os><ver> where 
# <os> is linux, mac or win and <ver> is the python version

# Note that on the Mac, by default the .so files are not portable because they 
# require libgcc_s.1.dylib and libgfortran.2.dylib to be installed at locations 
# defined by the compiler install. Three solutions are to patch 
# /Library/Frameworks/Python.framework/Versions/7.0/lib/python2.7/site-packages/numpy/distutils/fcompiler/gnu.py
#bht3:fcompiler bishop$ diff gnu.py gnu.py.original 
#128,129c128
#<             opt.extend(['-undefined', 'dynamic_lookup', '-bundle',
#<                         '-static-libgfortran','-static-libgcc'])
#---
#>             opt.extend(['-undefined', 'dynamic_lookup', '-bundle'])

# or (better) to define the env var
#  export LDFLAGS "-undefined dynamic_lookup -bundle -static-libgfortran -static-libgcc"

# or to edit path in the .so using:
# install_name_tool -change  <path>/libgcc_s.1.dylib @loader_path/libgcc_s.1.dylib <file>.so
# where is <path> is /usr/local/lib
# On Mac, use otool -L <file.so> to see what libraries are required

BIN = ../bin
LIBS = $(BIN)/pack_f.$(SUFFIX) $(BIN)/fellipse.$(SUFFIX) $(BIN)/pyspg.$(SUFFIX) $(BIN)/polymask.$(SUFFIX) $(BIN)/histogram2d.$(SUFFIX) $(BIN)/pytexture.$(SUFFIX)
SYMLIB := $(wildcard spsubs/*.for)
TXTLIB := $(wildcard texturesubs/*.for)
#----------------------------------------------------------------------
# linux (gfortran)
#COMPILER=--fcompiler=gnu95 
#PACKCOPTS=--f77flags="-fno-range-check"
#SUFFIX=so
#F2PY=f2py
#MOVE=mv
#DEL= \#
#----------------------------------------------------------------------
# mac (gfortran)
COMPILER=--fcompiler=gnu95 --f77exec=/sw/bin/gfortran
PACKCOPTS=--f77flags="-fno-range-check" # -static-libgfortran -static-libgcc
SUFFIX=so
F2PY=f2py
MOVE=mv
DEL=\#
RM=rm -f
#----------------------------------------------------------------------
# windows g77
#COMPILER=--fcompiler=gnu 
#PACKCOPTS=--f77flags="-fno-range-check"
#SUFFIX=pyd
#F2PY=f2py
#MOVE=copy
#DEL=del
#----------------------------------------------------------------------
# Windows gfortran-32
#COMPILER=--compiler=mingw32
#PACKCOPTS=--f77flags="-fno-range-check"
#SUFFIX=pyd
#F2PY=f2py
#MOVE=copy
#DEL=del




ask: 
	@echo ""
	@echo "Use make all or choose a target: "
	@echo "	$(LIBS) "

all:: $(BIN) $(LIBS)

# OSX: note that this is building .so's that require libgfortran and 
# libgcc_s.1 at runtime. Based on advice from Ilan at EPD, the .so files are 
# modified so that the libraries are placed in the same location as the .so files
#
$(BIN)::
	-mkdir $(BIN)

$(BIN)/pack_f.$(SUFFIX): pack_f.for
	$(F2PY) -c pack_f.for -m pack_f $(COMPILER) $(PACKCOPTS)
	$(MOVE) pack_f.$(SUFFIX) $(BIN)
	$(DEL) pack_f.$(SUFFIX)

$(BIN)/fellipse.$(SUFFIX): fellipse.for
	$(F2PY) -c fellipse.for -m fellipse $(COMPILER) $(PACKCOPTS)
	$(MOVE) fellipse.$(SUFFIX) $(BIN)
	$(DEL) fellipse.$(SUFFIX)

$(BIN)/polymask.$(SUFFIX): polymask.for
	$(F2PY) -c polymask.for -m polymask $(COMPILER)
	$(MOVE) polymask.$(SUFFIX) $(BIN)
	$(DEL) polymask.$(SUFFIX)

$(BIN)/histogram2d.$(SUFFIX): histogram2d.for
	$(F2PY) -c histogram2d.for -m histogram2d $(COMPILER)
	$(MOVE) histogram2d.$(SUFFIX) $(BIN)
	$(DEL) histogram2d.$(SUFFIX)

$(BIN)/pyspg.$(SUFFIX): pyspg.for $(SYMLIB)
	$(F2PY) -c pyspg.for $(SYMLIB) -m pyspg $(COMPILER)
	$(MOVE) pyspg.$(SUFFIX) $(BIN)
	$(DEL) pyspg.$(SUFFIX)

$(BIN)/pytexture.$(SUFFIX): $(TXTLIB)
	$(F2PY) -c $(TXTLIB) -m pytexture $(COMPILER)
	$(MOVE) pytexture.$(SUFFIX) $(BIN)
	$(DEL) pytexture.$(SUFFIX)


# basic outline for build; change ????? and add any specific options to F2PY step
# $(BIN)/?????.$(SUFFIX): ?????.for $(SYMLIB) $(BIN)
# 	$(F2PY) -c ?????.for -m ????? $(COMPILER)
# 	$(MOVE) ?????.$(SUFFIX) $(BIN)
# 	$(DEL) ?????.$(SUFFIX)

clean:
	$(RM) $(LIBS)
