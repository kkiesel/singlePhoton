#include <iostream>
#include <math.h>
#include <string>
#include <map>
#include <set>
#include <fstream>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TPRegexp.h"
#include "TArrayI.h"

#include "SusyEvent.h"
#include "TreeObjects.h"

class TreeWriter {
	public :
		TreeWriter(std::string inputName, std::string outputName, int loggingVerbosity_);
		TreeWriter(TChain* inputName, std::string outputName, int loggingVerbosity_);
		void Init( std::string outputName, int loggingVerbosity_ );
		virtual ~TreeWriter();
		virtual void Loop();
		bool isData();

		void SetProcessNEvents(int nEvents) { processNEvents = nEvents; }
		void SetReportEvents(unsigned int nEvents) { reportEvery = nEvents; }
		void SetLoggingVerbosity(unsigned int logVerb) { loggingVerbosity = logVerb; }
		void SetTriggerPaths( std::vector<const char*> const & tp ) { triggerNames = tp; }
		void PileUpWeightFile( std::string const & pileupFileName );
		void IncludeAJson( TString const & _fileName );

		TChain* inputTree;
		susy::Event* event;

		TFile* outFile;
		TTree* tree;
		TH1F* eventNumbers;

	private:
		bool passTrigger();
		bool isGoodLumi() const;
		float getPileUpWeight() const;

		int processNEvents; // number of events to be processed
		unsigned int reportEvery;
		unsigned int loggingVerbosity;
		// 0: no output
		// 1: only steps are shown
		// 2: object multiplicity shown
		// 3: detailed object info shown

		// important dataset information
		TH1F* pileupHisto;
		std::map<unsigned, std::set<unsigned> > goodLumiList;
		std::vector<const char*> triggerNames;


		// variables which will be stored in the tree
		std::vector<tree::Photon> photon;
		std::vector<tree::Jet> jet;
		std::vector<tree::Particle> electron;
		std::vector<tree::Particle> muon;
		std::vector<tree::Particle> genElectron;
		std::vector<tree::Particle> genPhoton;

		float met;
		float met_phi;
		float type1met;
		float type1met_phi;
		float ht;
		int nVertex;
		double weight;
		unsigned int runNumber;
		unsigned int eventNumber;
		unsigned int luminosityBlockNumber;
};
