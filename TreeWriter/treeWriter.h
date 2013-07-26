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
#include "TH2F.h"
#include "TH3I.h"
#include "TPRegexp.h"
#include "TArrayI.h"
#include "TLorentzVector.h"

#include "SusyEvent.h"
#include "TreeObjects.h"

class TreeWriter {
	public :
		TreeWriter(std::string inputName, std::string outputName, int loggingVerbosity_);
		TreeWriter(TChain* inputName, std::string outputName, int loggingVerbosity_);
		virtual ~TreeWriter();
		virtual void Loop();
		bool isData();

		void SetProcessNEvents(int nEvents) { processNEvents = nEvents; }
		void SetReportEvents(unsigned int nEvents) { reportEvery = nEvents; }
		void SetLoggingVerbosity(unsigned int logVerb) { loggingVerbosity = logVerb; }
		void SetTriggerPaths( std::vector<const char*> const & tp ) { triggerNames = tp; }
		void PileUpWeightFile( std::string const & pileupFileName );
		void IncludeAJson( TString const & _fileName );
		void UseAdditionalFOCut( bool v = true ) { useAdditionalFOCut = v; }

	private:
		void Init( std::string outputName, int loggingVerbosity_ );
		void SetBranches( TTree& tree );
		bool passTrigger();
		bool passRecommendedMetFilters() const;
		bool isGoodLumi() const;
		float getPileUpWeight() const;
		float getPtFromMatchedJet( const susy::Photon& myPhoton, bool isPhoton );
		float getHtHLT() const;
		float getSt( float ptCut, const std::vector<tree::Photon>& _photons ) const;
		std::vector<tree::Jet> getJets( const std::vector<tree::Photon>& excludedPhotons ) const;

		TChain* inputTree;
		susy::Event* event;

		// Saved objects
		TFile* outFile;
		TTree* photonTree;
		TTree* photonElectronTree;
		TTree* photonJetTree;
		TH1F* eventNumbers;
		std::map< std::string, TH2F* > hist2D;

		int processNEvents; // number of events to be processed
		unsigned int reportEvery;
		unsigned int loggingVerbosity;
		bool useAdditionalFOCut;

		// important dataset information
		TH1F* pileupHisto;
		std::map<unsigned, std::set<unsigned> > goodLumiList;
		std::vector<const char*> triggerNames;

		// variables which will be stored in the tree
		std::vector<tree::Photon> photons;
		std::vector<tree::Photon> photonElectrons;
		std::vector<tree::Photon> photonJets;
		std::vector<tree::Jet> jets;
		std::vector<tree::Particle> electrons;
		std::vector<tree::Particle> muons;
		std::vector<tree::Particle> genElectrons;
		std::vector<tree::Particle> genPhotons;

		float met;
		float type1met;
		float type0met;
		float htHLT;
		float st30;
		float st80;
		unsigned int nVertex;
		double weight;
		unsigned int runNumber;
		unsigned int eventNumber;
		unsigned int luminosityBlockNumber;
};
