void test(TString filename = ""){
	gSystem->Load("libSusyEvent.so");
	gSystem->Load("libTreeObjects.so");
	gROOT->LoadMacro("treeWriter.cc++");

	if (filename == ""){
		filename = "dcache:dcap://dcache-cms-dcap.desy.de/pnfs/desy.de/cms/tier2/store/user/kiesel/GJets_HT-400ToInf_8TeV-madgraph_v2/nTuplesSusy_GJets_V01/ef291a8fb6d60ba6695c261402fb808b/susyEvents_122_1_MKy.root";
		filename="../output/susyEvents_qcd.root";
	}

	TreeWriter *tw = new TreeWriter( filename, "susyTree.root", 0 );

	tw->PileUpWeightFile("pileUpReweighting/puWeights.root");

	tw->SetProcessNEvents(-1);
	tw->SkimEvents(true);
	tw->Loop();
}

