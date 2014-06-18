#include <iostream>
#include <sstream>
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

enum runTypes {
	kFullTree,
	kTree, // without systematic shifting, 33% computing time
	kGMSB,
	kSimplifiedModel
};

class TreeWriter {
	public :
		TreeWriter( int nFiles, char** fileList, std::string const& );
		virtual ~TreeWriter();
		virtual void Loop( int jecScale=0 );

		void SetRunType( runTypes t ) { runType = t; }

		// Command line output settings
		void SetProcessNEvents(int nEvents) { processNEvents = nEvents; }
		void SetReportEvents(unsigned int nEvents) { reportEvery = nEvents; }
		void SetLoggingVerbosity(unsigned int logVerb) { loggingVerbosity = logVerb; }

		// Set tigger and input Files
		void SetTriggerPaths( std::vector<const char*> const & tp ) { triggerNames = tp; }
		void SetJsonFile( TString const & filename );
		void SetQcdWeightHisto( TH2F histo ) { qcdWeightHisto = histo; }
		void SetPileUpWeightHisto( TH1F histo ) { pileupHisto = histo; }
		void SetPileUpWeightHistoUp( TH1F histo ) { pileupHistoUp = histo; }
		void SetPileUpWeightHistoDown( TH1F histo ) { pileupHistoDown = histo; }
		runTypes GetRunType() { return runType; }

	private:
		runTypes runType;
		void SetBranches( TTree& tree );
		bool passTrigger();
		bool isGoodLumi() const;
		eventType whichEventType( std::vector<tree::Photon>& photons,
			std::vector<tree::Photon>& photonElectrons,
			std::vector<tree::Photon>& photonJets );

		float getPileUpWeight();
		float getHt() const;
		TVector3 getMhtVector() const;
		TVector3 getRecoilVector( eventType eType ) const;
		void fillJets( int jecScale );
		void fillGenParticles();
		void fillLeptons();
		unsigned int countGoodJets();
		void getQcdWeights( float pt, float ht, float & qcdWeight, float & qcdWeightError );

		// Command line output settings
		unsigned int reportEvery;
		int processNEvents;
		unsigned int loggingVerbosity;

		// Additional information for producing the output
		TH1F pileupHisto;
		TH1F pileupHistoUp;
		TH1F pileupHistoDown;
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

		float metSig;
		float sumEt;

		float met;
		float metPhi;

		float metShiftxy;
		float metShiftxyPhi;

		float metWOzBoson;
		float metWOzBosonPhi;

		float met01corr;
		float met01corrPhi;

		float mht;
		float mhtPhi;

		float recoil;
		float recoilPhi;

		float ht;
		float weight;
		float weightPuUp;
		float weightPuDown;
		unsigned int nGoodJets;
		unsigned int nTracksPV;
		unsigned int nVertex;
		unsigned int runNumber;
		unsigned int eventNumber;
		unsigned int luminosityBlockNumber;
};
