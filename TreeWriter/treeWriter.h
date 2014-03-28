#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TPRegexp.h"

#include "SusyEvent.h"
#include "TreeObjects.h"

//#include "../../CMSSW/CMSSW_5_3_8_patch3/src/CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h" // to access the JEC scales

enum eventType {
	kPhotonEvent,
	kElectronEvent,
	kJetEvent
};

class TreeWriter {
	public :
		TreeWriter( int nFiles, char** fileList, std::string const& );
		virtual ~TreeWriter();
		virtual void Loop();

		// Command line output settings
		void SetProcessNEvents(int nEvents) { processNEvents = nEvents; }
		void SetReportEvents(unsigned int nEvents) { reportEvery = nEvents; }
		void SetLoggingVerbosity(unsigned int logVerb) { loggingVerbosity = logVerb; }

		// Configure output version
		void SplitTree( bool v = true ) { splitting = v; }
		void SkrinkTree( bool v = true ) { shrinkTree = v; }
		void FinalDistriputionsOnly( bool v = true ) { onlyMetPlots = v; }
		void ApplyHadronicSelection( bool v = true ) { hadronicSelection = v; }
		void SetPhotonPtThreshold( float th = 80 ) { photonPtThreshold = th; }

		// Set tigger and input Files
		void SetTriggerPaths( std::vector<const char*> const & tp ) { triggerNames = tp; }
		void SetPileUpWeightFile( std::string const & filename );
		void SetJsonFile( TString const & filename );
		void SetQcdWeightFile( std::string const & filename );

	private:
		void SetBranches( TTree& tree );
		bool passTrigger();
		bool isGoodLumi() const;
		eventType whichEventType( std::vector<tree::Photon>& photons,
			std::vector<tree::Photon>& photonElectrons,
			std::vector<tree::Photon>& photonJets );

		float getPileUpWeight();
		void getPtFromMatchedJet( tree::Photon& myPhoton, bool isPhoton, bool isPhotonJet, bool isPhotonElectron );
		void getPtFromMatchedJet1( tree::Photon& myPhoton, bool isPhoton, bool isPhotonJet, bool isPhotonElectron );
		float getHt() const;
		float getMht() const;
		void fillJets();
		unsigned int countGoodJets( bool clean );
		void getQcdWeights( float pt, float ht, float & qcdWeight, float & qcdWeightError );

		// Command line output settings
		unsigned int reportEvery;
		int processNEvents;
		unsigned int loggingVerbosity;

		// Configure output version
		bool splitting;
		bool onlyMetPlots;
		bool hadronicSelection;
		bool shrinkTree;
		float photonPtThreshold;

		// Additional information for producing the output
		TH1F pileupHisto;
		TH2F qcdWeightHisto;
		std::map<unsigned, std::set<unsigned> > goodLumiList;
		std::vector<const char*> triggerNames;

		TChain inputTree;
		susy::Event event;

		// Objects which can be saved to the file
		// photons: All tight photons (signal photons)
		// photonJets: All loose photons (qcd fake object)
		// photonElectrons: All photons with pixel seeds (ewk fake object)
		TFile outFile;
		TTree photonTree;
		TTree photonElectronTree;
		TTree photonJetTree;
		TH1F eventNumbers;
		TH3I nPhotons;
		std::map< std::string, TH2F > hist2D;
		std::map< std::string, TH1F > hist1D;

		// Variables which will be stored in the tree
		std::vector<tree::Photon> photons;
		std::vector<tree::Photon> photonElectrons;
		std::vector<tree::Photon> photonJets;
		std::vector<tree::Jet> jets;
		std::vector<tree::Particle> electrons;
		std::vector<tree::Particle> muons;
		std::vector<tree::Particle> genElectrons;
		std::vector<tree::Particle> genPhotons;

		float met;
		float metSig;
		float mht;
		float type1met;
		float type0met;
		float ht;
		float weight;
		unsigned int nGoodJets;
		unsigned int nTracksPV;
		unsigned int nVertex;
		unsigned int runNumber;
		unsigned int eventNumber;
		unsigned int luminosityBlockNumber;
};
