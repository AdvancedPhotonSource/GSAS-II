# build the compiled fortran codes needed by GSAS-II

BIN = ..\\bin
LIBS = $(BIN)/pack_f.$(SUFFIX) $(BIN)/pyspg.$(SUFFIX) $(BIN)/polymask.$(SUFFIX) $(BIN)/histogram2d.$(SUFFIX) $(BIN)/pypowder.$(SUFFIX)
SYMLIB := $(wildcard spsubs/*.for)
POWLIB := $(wildcard powsubs/*.for)
MACFIX = #
MACCOPY = #
#----------------------------------------------------------------------
# linux (gfortran)
#GSASlib = ./libgsas.a
#COMPILER=--fcompiler=gnu95 
#PACKCOPTS=--f77flags="-fno-range-check"
#SUFFIX=so
#F2PY=f2py
#MOVE=mv
#DEL= \#
#----------------------------------------------------------------------
# mac (gfortran)
#GSASlib = /Users/toby/software/work/gsas/2009Aug31/libgsas.a
#COMPILER=--fcompiler=gnu95 --f90exec=/usr/local/bin/gfortran
#PACKCOPTS=--f77flags="-fno-range-check"
#SUFFIX=so
#F2PY=f2py
#MOVE=mv
#DEL=\#
#MACFIX = install_name_tool -change 
#MACCOPY = cp
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
COMPILER=--compiler=mingw32
PACKCOPTS=--f77flags="-fno-range-check"
SUFFIX=pyd
F2PY=f2py
MOVE=copy
DEL=del




ask: 
	@echo ""
	@echo "Use make all or choose a target: "
	@echo "	$(LIBS) "

all:: $(BIN) $(LIBS)

# OSX: note that this is building .so's that require libgfortran and 
# libgcc_s.1 at runtime. Based on advice from Ilan at EPD, the .so files are 
# modified so that the libraries are placed in the same location as the .so files
# Use otool -L <file.so> to see what is required
#
$(BIN)::
	-mkdir $(BIN)

$(BIN)/pack_f.$(SUFFIX): pack_f.for $(BIN)
	$(F2PY) -c pack_f.for -m pack_f $(COMPILER) $(PACKCOPTS)
#	$(MACFIX) /usr/local/lib/libgcc_s.1.dylib @loader_path/libgcc_s.1.dylib pack_f.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgcc_s.1.dylib $(BIN)
#	$(MACFIX) /usr/local/lib/libgfortran.2.dylib @loader_path/libgfortran.2.dylib pack_f.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgfortran.2.dylib $(BIN)
	$(MOVE) pack_f.$(SUFFIX) $(BIN)
	$(DEL) pack_f.$(SUFFIX)

$(BIN)/pypowder.$(SUFFIX): pypowder.for $(POWLIB) $(BIN)
	$(F2PY) -c pypowder.for -m pypowder $(COMPILER) $(POWLIB)
#	$(MACFIX) /usr/local/lib/libgcc_s.1.dylib @loader_path/libgcc_s.1.dylib pypowder.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgcc_s.1.dylib $(BIN)
#	$(MACFIX) /usr/local/lib/libgfortran.2.dylib @loader_path/libgfortran.2.dylib pypowder.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgfortran.2.dylib $(BIN)
	$(MOVE) pypowder.$(SUFFIX) $(BIN)
	$(DEL) pypowder.$(SUFFIX)

$(BIN)/pyspg.$(SUFFIX): pyspg.for $(SYMLIB) $(BIN)
	$(F2PY) -c pyspg.for $(SYMLIB) -m pyspg $(COMPILER)
#	$(MACFIX) /usr/local/lib/libgcc_s.1.dylib @loader_path/libgcc_s.1.dylib pyspg.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgcc_s.1.dylib $(BIN)
#	$(MACFIX) /usr/local/lib/libgfortran.2.dylib @loader_path/libgfortran.2.dylib pyspg.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgfortran.2.dylib $(BIN)
	$(MOVE) pyspg.$(SUFFIX) $(BIN)
	$(DEL) pyspg.$(SUFFIX)

$(BIN)/polymask.$(SUFFIX): polymask.for $(SYMLIB) $(BIN)
	$(F2PY) -c polymask.for -m polymask $(COMPILER)
#	$(MACFIX) /usr/local/lib/libgcc_s.1.dylib @loader_path/libgcc_s.1.dylib polymask.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgcc_s.1.dylib $(BIN)
#	$(MACFIX) /usr/local/lib/libgfortran.2.dylib @loader_path/libgfortran.2.dylib polymask.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgfortran.2.dylib $(BIN)
	$(MOVE) polymask.$(SUFFIX) $(BIN)
	$(DEL) polymask.$(SUFFIX)

$(BIN)/histogram2d.$(SUFFIX): histogram2d.for $(SYMLIB) $(BIN)
	$(F2PY) -c histogram2d.for -m histogram2d $(COMPILER)
#	$(MACFIX) /usr/local/lib/libgcc_s.1.dylib @loader_path/libgcc_s.1.dylib histogram2d.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgcc_s.1.dylib $(BIN)
#	$(MACFIX) /usr/local/lib/libgfortran.2.dylib @loader_path/libgfortran.2.dylib histogram2d.$(SUFFIX)
#	$(MACCOPY) /usr/local/lib/libgfortran.2.dylib $(BIN)
	$(MOVE) histogram2d.$(SUFFIX) $(BIN)
	$(DEL) histogram2d.$(SUFFIX)


# basic outline for build; change ????? and add any specific options to F2PY step
# $(BIN)/?????.$(SUFFIX): ?????.for $(SYMLIB) $(BIN)
# 	$(F2PY) -c ?????.for -m ????? $(COMPILER)
# 	$(MACFIX) /usr/local/lib/libgcc_s.1.dylib @loader_path/libgcc_s.1.dylib ?????.$(SUFFIX)
# 	$(MACCOPY) /usr/local/lib/libgcc_s.1.dylib $(BIN)
# 	$(MACFIX) /usr/local/lib/libgfortran.2.dylib @loader_path/libgfortran.2.dylib ?????.$(SUFFIX)
# 	$(MACCOPY) /usr/local/lib/libgfortran.2.dylib $(BIN)
# 	$(MOVE) ?????.$(SUFFIX) $(BIN)
# 	$(DEL) ?????.$(SUFFIX)


# no longer in use
#$(BIN)/fitellipse.$(SUFFIX): fitellipse.for
#	cd $(BIN); $(F2PY) -c ../fitellipse.for -m fitellipse --fcompiler=gfortran --f90exec=/usr/local/bin/gfortran --f77flags="-fno-range-check"


