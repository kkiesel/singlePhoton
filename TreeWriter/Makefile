LIBS = $(shell root-config --libs)
INCS = -I$(shell root-config --incdir) -I.

WARN = -Wall -Wshadow

SHARED_LIB = libSusyEvent.so libTreeObjects.so
OBJ = treeWriter.o main.o
EXE = treeWriter

all: $(SHARED_LIB) $(EXE) Makefile
	@if [ -d /afs ]; then make -C pileUpReweighting --no-print-directory; fi

%.o:%.cc
	g++ -std=gnu++11 -o $@ $+ -c -O2 $(INCS) $(WARN)

$(EXE): $(OBJ) $(SHARED_LIB)
	g++ -std=gnu++11 -O2 -o $@ $+ $(LIBS) $(WARN)


# see http://root.cern.ch/drupal/content/interacting-shared-libraries-rootcint
# The libaries will be copied to your cmssw lib directory. If this is not possible,
# please execute 'LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`'
lib%.so: %.h %_LinkDef.h
	rootcint -f $(subst so,cxx,$@) -c $+
	g++ -std=gnu++11 -O2 -Wall  -shared -fPIC `root-config --ldflags` $(INCS) -o $@ $(subst lib,,$(subst so,cc,$@)) $(subst so,cxx,$@)
	@if [ -d $(CMSSW_BASE)/lib/$(SCRAM_ARCH) -a "$(CMSSW_BASE)" != "" ]; then\
		cp $@ $(CMSSW_BASE)/lib/$(SCRAM_ARCH);\
	else\
		echo -e "\n#############################################\n### Please copy the library manually to $(root-config --libdir)! ### \n#############################################\n";\
	fi

clean:
	@rm -f $(OBJ) $(SHARED_LIB)
	@rm -f *.cxx lib*h # remove rootcints generated source and header files
	@rm -f AutoDict_* *_h.d *_cc.d *_C.d # cint generated files


