LIBS = $(shell root-config --libs)
INCS = -I$(shell root-config --incdir) -I.

TARGETS = libSusyEvent.so libTreeObjects.so

all: $(TARGETS)
	make -C pileUpReweighting

# see http://root.cern.ch/drupal/content/interacting-shared-libraries-rootcint
lib%.so: %.h %_LinkDef.h
	rootcint -f $(subst so,cxx,$@) -c $+
	g++ -O2 -Wall  -shared -fPIC $(INCS) $(LIBS) -o $@ $(subst lib,,$(subst so,cc,$@)) $(subst so,cxx,$@)

clean:
	@rm -f *.o *.so
	@rm -f *.cxx lib*h # remove rootcints generated source and header files
	@rm -f AutoDict_* *_h.d *_cc.d *_C.d # cint generated files

