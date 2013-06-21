#include "treeWriter.h"

char* getCmdOption(char ** begin, char ** end, const std::string & option) {
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
		return *itr;
	return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
	return std::find(begin, end, option) != end;
}

void printHelp() {
	std::cout << "Usage: ./executable [-h] [-MC] -o outFileName.root -i inputFileName1.root inputFileName2.root ..." << std::endl;
}

int main( int argc, char** argv ) {

	if( cmdOptionExists( argv, argv+argc, "-h" ) ) {
		printHelp();
		return 0;
	}

	bool isMC = cmdOptionExists( argv, argv+argc, "-MC" ); // false by default

	char * filename = getCmdOption(argv, argv + argc, "-o");
	std::string outputFileName;
	if (filename) {
		outputFileName = filename;
	} else {
		std::cout << "No output filename specified." << std::endl;
		printHelp();
		return 1;
	}

	TChain *inputTree = new TChain("susyTree");
	char ** itr = std::find(argv, argv+argc, (std::string)"-i");
	while (itr != argv && ++itr != argv+argc)
		inputTree->Add( *itr );

	std::cout << "Write to output file \"" << outputFileName << "\"" << std::endl;

	TreeWriter *tw = new TreeWriter( inputTree, outputFileName, 0 );
	// common settings
	tw->SetProcessNEvents(-1);
	tw->SetReportEvents(20000);

	if( isMC ) {
		tw->PileUpWeightFile("pileUpReweighting/puWeights.root");
	} else {
		tw->IncludeAJson("../../Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt");
		std::vector<const char*> triggerNames;
		triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFHT400_v" );
		triggerNames.push_back( "HLT_Photon70_CaloIdXL_PFNoPUHT400_v" );
		tw->SetTriggerPaths( triggerNames );
	}

	double start_time = time(NULL);
	tw->Loop();
	double end_time = time(NULL);

	std::cout << "Job needed " << 1.*(end_time - start_time)/3600 << " h real time." << std::endl;
	return 0;
}

