#include "treeWriter.h"

int main( int argc, char** argv ) {

	if( argc < 2 ) {
		std::cout << "usage: ./execute outputFileName.root input1.root input2.root ..." << std::endl;
		return 1;
	}
	std::string outputFileName = argv[1];

	TChain *inputTree = new TChain("susyTree");
	for( int i=2; i<argc; ++i)
		inputTree->Add( argv[i] );
	TreeWriter *tw = new TreeWriter( inputTree, outputFileName, 0 );

	// common settings
	tw->SplitTree( true );
	//tw->FinalDistriputionsOnly();
	tw->SetPhotonPtThreshold( 80 );
	tw->ApplyHadronicSelection( true );

	tw->SetQcdWeightFile("../plotTree/qcdWeight.root");
	tw->SetJsonFile( "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" );
	tw->SetJsonFile( "../../Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" );
	tw->SetPileUpWeightFile( "pileUpReweighting/puWeights.root" );

	std::vector<const char*> triggerNames;
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFHT400_v" );
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFNoPUHT400_v" );
	tw->SetTriggerPaths( triggerNames );

	double start_time = time(NULL);
	tw->Loop();
	double end_time = time(NULL);

	std::cout << "Job needed " << 1.*(end_time - start_time)/3600 << " h real time." << std::endl;
	return 0;
}

