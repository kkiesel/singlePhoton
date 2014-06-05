#include "treeWriter.h"

template <class HIST>
HIST getHisto( std::string const & filename, std::string const & histname ) {
	TFile file( filename.c_str() );
	if( file.IsZombie() )
		std::cerr << "ERROR: Could not read pileup weight file " << filename << std::endl;

	if( file.GetListOfKeys()->Contains( histname.c_str() ) )
		return *((HIST*) file.Get( histname.c_str() ));
	else
		std::cerr << "ERROR: Could not extract " << histname
			<< " histogram from " << filename << std::endl;
	return HIST();
}

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
	std::string pileupScenario = "S10";

	// for fake rate studies
	//tw.SetPhotonPtThreshold( 0 ); tw.ApplyHadronicSelection( false ); tw.SkrinkTree( true );

	const std::string lumiJsonName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
	if( access( lumiJsonName.c_str(), F_OK ) != -1 )
		tw.SetJsonFile( lumiJsonName );
	else
		tw.SetJsonFile( "../../Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" );

	gSystem->Load("libHistPainter"); // to avoid waring and errors when reading th2 from file
	tw.SetQcdWeightHisto( getHisto<TH2F>( "../plotTree/qcdWeight.root", "qcdWeight" ) );
	tw.SetPileUpWeightHisto( getHisto<TH1F>( "pileUpReweighting/puWeights.root", "pileupWeight"+pileupScenario ) );
	tw.SetPileUpWeightHistoUp( getHisto<TH1F>( "pileUpReweighting/puWeights.root", "pileupWeightUp"+pileupScenario ) );
	tw.SetPileUpWeightHistoDown( getHisto<TH1F>( "pileUpReweighting/puWeights.root", "pileupWeightDown"+pileupScenario ) );

	std::vector<const char*> triggerNames;
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFHT400_v" );
	triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFNoPUHT400_v" );
	tw.SetTriggerPaths( triggerNames );

	double start_time = time(NULL);
	tw.Loop();
	tw.Loop(1);
	tw.Loop(-1);
	double end_time = time(NULL);

	std::cout << "Job needed " << 1.*(end_time - start_time)/3600 << " h real time." << std::endl;
	return 0;
}

