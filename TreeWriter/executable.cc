#include "treeWriter.h"

int main( int argc, char** argv ) {

	if( argc < 3 ) {
		std::cout << "usage: ./execute outputFileName.root input1.root input2.root ..." << std::endl;
		return 1;
	}

	TreeWriter tw( argc-2, argv+2, argv[1] );

	// common settings
	tw.SplitTree( true );
	tw.FinalDistriputionsOnly();
	tw.SetPhotonPtThreshold( 110 );
	tw.ApplyHadronicSelection( true );

	tw.SetQcdWeightFile("qcdWeight.root");
	tw.SetJsonFile( "Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" );
	tw.SetPileUpWeightFile( "puWeights.root" );
	tw.SetPileUpWeightFile( "pileUpReweighting/puWeights.root" );

	std::vector<const char*> triggerNames;
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFHT400_v" );
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFNoPUHT400_v" );
	tw.SetTriggerPaths( triggerNames );

	double start_time = time(NULL);
	tw.Loop();
	double end_time = time(NULL);

	std::cout << "Job needed " << 1.*(end_time - start_time)/3600 << " h real time." << std::endl;
	return 0;
}

