#include<iostream>
#include<math.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"

#include "SusyEvent.h"
#include "TreeObjects.h"

class TreeWriter {
	public :
		TreeWriter(TString inputName, TString outputName, int loggingVerbosity_);
		TreeWriter(TChain* inputName, TString outputName, int loggingVerbosity_);
		void Init( TString outputName, int loggingVerbosity_ );
		virtual ~TreeWriter();
		virtual void Loop();

		void SetProcessNEvents(int nEvents) { processNEvents = nEvents; }
		void SetReportEvents(int nEvents) { reportEvery = nEvents; }
		void SetLoggingVerbosity(int logVerb) { loggingVerbosity = logVerb; }
		void SkimEvents(bool skim_){ skim = skim_; }
		void PileUpWeightFile( string pileupFileName );

		TChain *inputTree;
		susy::Event *event;

		TFile *outFile;
		TTree *tree;
		TH1F *eventNumbers;

		float getPtFromMatchedJet( susy::Photon, susy::Event );
		float deltaR( TLorentzVector, TLorentzVector );

	private:
		int processNEvents; // number of events to be processed
		int reportEvery;
		int loggingVerbosity;
		bool skim;

		// important dataset information
		TH1F* pileupHisto;

		// variables which will be stored in the tree
		std::vector<tree::Photon> photon;
		std::vector<tree::Jet> jet;
		std::vector<tree::Particle> electron;
		std::vector<tree::Particle> muon;

		float met;
		float met_phi;
		float type1met;
		float type1met_phi;
		float ht;
		int nVertex;
		float pu_weight;
};

