#include "treeWriter.h"

int main( int argc, char** argv ) {

	std::string filename = "../susyEvents.root";
	std::string outputFilename = "susyTree.root";

	TreeWriter *tw = new TreeWriter( filename, outputFilename, 0 );

	// settings
	tw->PileUpWeightFile("pileUpReweighting/puWeights.root");
	tw->SetProcessNEvents(10);
	tw->SkimEvents(true);
	tw->Loop();

}

