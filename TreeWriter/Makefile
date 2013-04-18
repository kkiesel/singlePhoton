LIBS = $(shell root-config --libs)
INCS = -I$(shell root-config --incdir) -I.
#LIBS = /cvmfs/cms.cern.ch/slc5_amd64_gcc462/lcg/root/5.34.02-cms/lib
#INCS = -I/cvmfs/cms.cern.ch/slc5_amd64_gcc462/lcg/root/5.34.02-cms/include

SHARED_LIB = libSusyEvent.so libTreeObjects.so
OBJ = treeWriter.o executable.o
EXE = executable

all: $(SHARED_LIB) $(EXE)
	@make -C pileUpReweighting --no-print-directory

%.o:%.cc
	g++ -o $@ $+ -c $(INCS)

$(EXE): $(OBJ) $(SHARED_LIB)
	g++ -o $@ $+ $(LIBS) && ./$@


# see http://root.cern.ch/drupal/content/interacting-shared-libraries-rootcint
lib%.so: %.h %_LinkDef.h
	rootcint -f $(subst so,cxx,$@) -c $+
	g++ -O2 -Wall  -shared -fPIC `root-config --ldflags` $(INCS) -o $@ $(subst lib,,$(subst so,cc,$@)) $(subst so,cxx,$@)
	cp $@ $(CMSSW_BASE)/lib/$(SCRAM_ARCH)

clean:
	@rm -f $(OBJ) $(SHARED_LIB)
	@rm -f *.cxx lib*h # remove rootcints generated source and header files
	@rm -f AutoDict_* *_h.d *_cc.d *_C.d # cint generated files


