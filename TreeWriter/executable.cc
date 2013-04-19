#include "treeWriter.h"

int main( int argc, char** argv ) {

	if( argc < 2 ) {
		std::cout << "usage: ./execute outputFileName.root input1.root input2.root ..." << std::endl;
		return 1;
	}
	std::string outputFileName = argv[1];

	TChain *inputTree = new TChain("susyTree");
	for( unsigned int i=2; i<argc; ++i)
		inputTree->Add( argv[i] );

	std::cout << "Write to output file \"" << outputFileName << "\"" << std::endl;

	TreeWriter *tw = new TreeWriter( inputTree, outputFileName, 0 );

	// settings
	tw->PileUpWeightFile("pileUpReweighting/puWeights.root");
	tw->SetProcessNEvents(-1);
	tw->SetReportEvents(1000);
	tw->SkimEvents(true);
	tw->Loop();

}

