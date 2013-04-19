#include "treeWriter.h"

int main( int argc, char** argv ) {

	std::string filename = "../susyEvents.root";
	//std::string filename = "dcache:dcap://dcache-cms-dcap.desy.de/pnfs/desy.de/cms/tier2/store/user/jschulz/nTuples/WJets_V01/susyEvents_812_1_Fbi.root";
	std::string outputFilename = "susyTree.root";

	TreeWriter *tw = new TreeWriter( filename, outputFilename, 0 );

	// settings
	tw->PileUpWeightFile("pileUpReweighting/puWeights.root");
	tw->SetProcessNEvents(-1);
	tw->SetReportEvents(1000);
	tw->SkimEvents(true);
	tw->Loop();

}

