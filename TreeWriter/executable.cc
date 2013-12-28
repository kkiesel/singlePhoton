#include "treeWriter.h"

int main( int argc, char** argv ) {

	if( argc < 3 ) {
		std::cout << "usage: ./execute outputFileName.root input1.root input2.root ..." << std::endl;
		return 1;
	}

	TreeWriter tw( argc-2, argv+2, argv[1] );

	// common settings
	tw.SplitTree( true );
	tw.FinalDistriputionsOnly( false );
	tw.SetPhotonPtThreshold( 110 );
	tw.ApplyHadronicSelection( true );

	// for fake rate studies
	tw.SetPhotonPtThreshold( 0 );
	tw.ApplyHadronicSelection( false );
	tw.SkrinkTree( true );

	tw.SetQcdWeightFile("../plotTree/qcdWeight.root");
	const std::string lumiJsonName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
	if( access( lumiJsonName.c_str(), F_OK ) != -1 )
		tw.SetJsonFile( lumiJsonName );
	else
		tw.SetJsonFile( "../../Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" );
	tw.SetPileUpWeightFile( "pileUpReweighting/puWeights.root" );

	std::vector<const char*> triggerNames;
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFHT400_v" );
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFNoPUHT400_v" );
	tw.SetTriggerPaths( triggerNames );

	double start_time = time(NULL);
	tw.Loop();
	//tw.LoopElectronsOnly();
	double end_time = time(NULL);

	std::cout << "Job needed " << 1.*(end_time - start_time)/3600 << " h real time." << std::endl;
	return 0;
}

