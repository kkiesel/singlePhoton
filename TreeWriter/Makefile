# convention: name all shared cint objects with .so, all other objects with .o

LIBS = $(shell root-config --libs)
INCS = -I$(shell root-config --incdir) -I.

OBJ = treeWriter.o TreeObjects.o SusyEvent.o EventDict.so TreeDict.so

treeWriter: $(OBJ)
	g++ $+ -o $@ $(LIBS)

%.so: %.cxx
	g++ $(INCS) -c -shared -fPIC -o $@ $<

%.o: %.cc
	g++ $(INCS) -c -o $@ $<

EventDict.cxx: SusyEvent.h SusyNtuplizer_LinkDef.h
	rootcint -f $@ -c $+

TreeDict.cxx: TreeObjects.h TreeObjects_LinkDef.h
	rootcint -f $@ -c $+


clean:
	@rm -f $(subst .so,.h,$(OBJ)) *.so *.cxx treeWriter
