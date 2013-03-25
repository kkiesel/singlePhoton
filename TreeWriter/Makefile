LIBS = $(shell root-config --libs)
INCS = -I$(shell root-config --incdir) -I.

OBJ = treeWriter.o SusyEvent.o TreeObjects.o libSusyEvent.so libTreeObjects.so

#TODO: set library path correctly, eg in your .bashrc
treeWriter: $(OBJ)
	g++ $+ -o $@ $(LIBS)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
	./$@

%.o: %.cc
	g++ $(INCS) -c -o $@ $<

# see http://root.cern.ch/drupal/content/interacting-shared-libraries-rootcint
lib%.so: %.h %_LinkDef.h
	rootcint -f $(subst so,cxx,$@) -c -p $+
	g++ $(INCS) -shared -fPIC -o $@ $(subst so,cxx,$@) $(subst lib,,$(subst so,cc,$@))

clean:
	@rm -f *.o *.so
	@rm -f *.cxx lib*h # remove rootcints generated source and header files
	@rm -f treeWriter # remove executable

