
// syntax: root -b -q -l "test.C(\"input.root\",\"output.root\")"

void test(string filename="", string outputFilename=""){

	// set default parameters
	if( filename == "" )
		filename = "dcache:dcap://dcache-cms-dcap.desy.de/pnfs/desy.de/cms/tier2/store/user/kiesel/GJets_HT-400ToInf_8TeV-madgraph_v2/nTuplesSusy_GJets_V01/ef291a8fb6d60ba6695c261402fb808b/susyEvents_122_1_MKy.root";
	if( outputFilename == "" )
		outputFilename = "susyTree.root";

	// if input file is no a root file, a directory with root files is assumed
	if( filename.size() - filename.rfind(".root") - 5 ) {
		cout << "this is a directory" << endl;
	}

	// load stuff
	gSystem->Load("libSusyEvent.so");
	gSystem->Load("libTreeObjects.so");
	gROOT->LoadMacro("treeWriter.cc++");
	TreeWriter *tw = new TreeWriter( filename, outputFilename, 0 );

	// settings
	tw->PileUpWeightFile("pileUpReweighting/puWeights.root");
	tw->SetProcessNEvents(5);
	tw->SkimEvents(true);
	tw->Loop();

}

