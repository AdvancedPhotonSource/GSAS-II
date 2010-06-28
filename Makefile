# build the compiled fortran codes needed by GSAS-II

BIN = bin
LIBS = $(BIN)/pack_f.$(SUFFIX) $(BIN)/pyspg.$(SUFFIX) 
LIBSwGSAS = $(BIN)/pypowder.$(SUFFIX)
SYMLIB := $(wildcard spsubs/*.for)
#----------------------------------------------------------------------
# linux (gfortran)
#COMPILER=--fcompiler=gnu95 
#PACKCOPTS=--f77flags="-fno-range-check"
#SUFFIX=so
#F2PY=f2py
#MOVE=mv
#DEL=echo
#----------------------------------------------------------------------
# mac (gfortran)
GSASlib = /Users/toby/software/work/gsas/2009Aug31/libgsas.a
COMPILER=--fcompiler=gnu95 --f90exec=/usr/local/bin/gfortran
#PACKCOPTS=--f77flags="-fno-range-check -static-libgcc"
SUFFIX=so
F2PY=f2py
MOVE=mv
DEL=echo
#----------------------------------------------------------------------
# windows g77
#COMPILER=--fcompiler=gnu 
#PACKCOPTS=--f77flags="-fno-range-check"
#SUFFIX=pyd
#F2PY=f2py.py
#MOVE=copy
#DEL=del
#----------------------------------------------------------------------

ask: 
	@echo ""
	@echo "Use make all or choose a target: "
	@echo "	$(LIBS) $(LIBSwGSAS)"
	@echo "   Note: target $(LIBSwGSAS) requires the GSAS object library."
	@echo "     This is not built with make all. You will need to edit the"
	@echo "     Makefile to set GSASlib to point to the correct location."

all:: $(BIN) $(LIBS)

# OSX: note that this is building .so's that require libgfortran and 
# libgcc_s.1 at runtime. The former can be removed by renaming the accessed 
# libgfortran.3.dynlib. Not sure how to avoid libgcc_s.1
# Use otool -L <file.so> to see what is required
#
.PHONY: $(BIN)
	mkdir $(BIN)

$(BIN)/pack_f.$(SUFFIX): pack_f.for
	$(F2PY) -c pack_f.for -m pack_f $(COMPILER) $(PACKCOPTS)
	$(MOVE) pack_f.$(SUFFIX) $(BIN)
	$(DEL) pack_f.$(SUFFIX)

$(BIN)/pypowder.$(SUFFIX): pypowder.for $(GSASlib)
	$(F2PY) -c pypowder.for -m pypowder $(COMPILER) $(GSASlib)
	$(MOVE) pypowder.$(SUFFIX) $(BIN)
	$(DEL) pypowder.$(SUFFIX)

$(BIN)/pyspg.$(SUFFIX): pyspg.for $(SYMLIB)
	$(F2PY) -c pyspg.for $(SYMLIB) -m pyspg $(COMPILER) 
	$(MOVE) pyspg.$(SUFFIX) $(BIN)
	$(DEL) pyspg.$(SUFFIX)

# no longer in use
#$(BIN)/fitellipse.$(SUFFIX): fitellipse.for
#	cd $(BIN); $(F2PY) -c ../fitellipse.for -m fitellipse --fcompiler=gfortran --f90exec=/usr/local/bin/gfortran --f77flags="-fno-range-check"


