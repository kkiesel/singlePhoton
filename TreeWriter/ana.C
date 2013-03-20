void ana(){
	gSystem->Load("libSusyEvent.so");
	gROOT->LoadMacro("treeWriter.cc+");

	TreeWriter *tw = new TreeWriter("qcd-1000-nTuple-test.root", "myTree.root");
	tw->SetProcessNEvents(10);
	tw->Loop();
}
