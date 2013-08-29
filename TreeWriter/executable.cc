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

	std::cout << "Write to output file \"" << outputFileName << "\"" << std::endl;

	TreeWriter *tw = new TreeWriter( inputTree, outputFileName, 0 );
	bool isData = tw->isData();

	// common settings
	tw->SplitTree( false );
	tw->SetPhotonPtThreshold( 80 );
	tw->ApplyHadronicSelection( false );

	if( isData ) {
		std::cout << "Process data." << std::endl;

		// try json of afs first, else try out reading local one
		if( ! tw->IncludeAJson( "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" )
				&& ! tw->IncludeAJson( "../../Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" ) ) {
			std::cout << "Can't read either of the two json files, abort." << std::endl;
			return 1;
		}
		std::vector<const char*> triggerNames;
		triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFHT400_v" );
		triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFNoPUHT400_v" );
		tw->SetTriggerPaths( triggerNames );
	} else {
		std::cout << "Process simulation." << std::endl;
		tw->PileUpWeightFile("pileUpReweighting/puWeights.root");
	}

	double start_time = time(NULL);
	tw->Loop();
	double end_time = time(NULL);

	std::cout << "Job needed " << 1.*(end_time - start_time)/3600 << " h real time." << std::endl;
	return 0;
}

