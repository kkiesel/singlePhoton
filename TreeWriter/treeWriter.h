#include<iostream>
#include<math.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"

#include "SusyEvent.h"
#include "TreeObjects.h"

class TreeWriter {
	public :
		TreeWriter(TString inputName, TString outputName, int loggingVerbosity_);
		virtual ~TreeWriter();
		virtual void Loop();

		void SetProcessNEvents(int nEvents) { processNEvents = nEvents; }
		void SetReportEvents(int nEvents) { reportEvery = nEvents; }
		void SetLoggingVerbosity(int logVerb) { loggingVerbosity = logVerb; }
		void SkimEvents(bool skim_){ skim = skim_; }

		TChain *inputTree;
		susy::Event *event;

		TFile *outFile;
		TTree *tree;

		float getPtFromMatchedJet( susy::Photon, susy::Event );
		float deltaR( TLorentzVector, TLorentzVector );

	private:
		int processNEvents; // number of events to be processed
		int reportEvery;
		int loggingVerbosity;
		bool skim;

		// variables which will be stored in the tree
		std::vector<tree::Photon> photon;
		std::vector<tree::Jet> jet;
		float met;
		float ht;
		int nVertex;
		int nElectron;
};

