void test(){
	gSystem->Load("libSusyEvent.so");
	gSystem->Load("libTreeObjects.so");
	gROOT->LoadMacro("treeWriter.cc++");

	TString filename = "dcache:dcap://dcache-cms-dcap.desy.de/pnfs/desy.de/cms/tier2/store/user/kiesel/GJets_HT-400ToInf_8TeV-madgraph_v2/nTuplesSusy_GJets_V01/ef291a8fb6d60ba6695c261402fb808b/susyEvents_122_1_MKy.root";

	TreeWriter *tw = new TreeWriter( filename, "myTree.root", 3);
	tw->SetProcessNEvents(1);
	tw->SkimEvents(true);
	tw->SetLoggingVerbosity(3);
	tw->Loop();
}

