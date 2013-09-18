#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TPRegexp.h"

#include "SusyEvent.h"
#include "TreeObjects.h"

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
		float getPileUpWeight();
		float getPtFromMatchedJet( susy::Photon const & myPhoton, bool isPhoton );
		float getHtHLT() const;
		float getHt() const;
		float getJetHt() const;
		std::vector<tree::Jet> getJets( bool clean ) const;
		void getQcdWeights( float pt, float ht, float & qcdWeight, float & qcdWeightUp, float & qcdWeightDown );

		// Command line output settings
		unsigned int reportEvery;
		int processNEvents;
		unsigned int loggingVerbosity;

		// Configure output version
		bool splitting;
		bool onlyMetPlots;
		bool hadronicSelection;
		float photonPtThreshold;

		// Additional information for producing the output
		TH1F pileupHisto;
		TH2F qcdWeightHisto;
		std::map<unsigned, std::set<unsigned> > goodLumiList;
		std::vector<const char*> triggerNames;

		TChain inputTree;
		susy::Event event;

		// Objects which can be saved to the file
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
		float type1met;
		float type0met;
		float htHLT;
		float ht;
		float jetHt;
		float weight;
		unsigned int nVertex;
		unsigned int runNumber;
		unsigned int eventNumber;
		unsigned int luminosityBlockNumber;
};
