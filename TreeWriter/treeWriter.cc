#include<iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "SusyEvent.h"


namespace tree {
// In this namespace classes for the trees are defined.

class Particle {
	public:
		float pt, eta, phi;
};

class Photon : public Particle {
	public:
		float r9, sigmaIetaIeta, hadTowOverEm, pixelseed;
		float chargedIso, neutralIso, photonIso;
};

class Jet : public Particle{
	public:
		float pt, eta;
};

bool EtGreater(const tree::Particle p1, const tree::Particle p2) {
  return p1.pt > p2.pt;
}

} // end namespace definition

class TreeWriter {
	public :
		TreeWriter(TString inputName, TString outputName );
		virtual ~TreeWriter();
		virtual void Loop();

		void SetProcessNEvents(int nEvents) { processNEvents = nEvents; }
		void SetReportEvents(int nEvents) { reportEvery = nEvents; }

		TFile *inputFile;
		TTree *inputTree;
		susy::Event *event;

		TFile *outFile;
		TTree *tree;

	private:
		int processNEvents; // number of events to be processed
		int reportEvery;
		int loggingVerbosity;

		// variables which will be stored in the tree
		std::vector<tree::Photon> photon;
		std::vector<tree::Jet> jet;
		float met;
		int nVertex;
		float weight;
};


TreeWriter::TreeWriter(TString inputName, TString outputName) {
	// read the input file
	inputFile = new TFile( inputName, "read" );
	inputTree = (TTree*) inputFile->Get("susyTree");
	event = new susy::Event;
	inputTree->SetBranchAddress("susyEvent", &event);

	// open the output file
	outFile = new TFile( outputName, "recreate" );
	tree = new TTree("susyTree","Tree for single photon analysis");

	// set default parameter
	processNEvents = -1;
	reportEvery = 1000;
	loggingVerbosity = 0;

}

TreeWriter::~TreeWriter() {
	if (!inputTree) return;
	delete inputTree->GetCurrentFile();
}

void TreeWriter::Loop() {
	// here the event loop is implemented and the tree is filled
	if (inputTree == 0) return;

	// get number of events to be proceeded
	Long64_t nentries = inputTree->GetEntries();
	if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

	if( loggingVerbosity > 0 )
		std::cout << "Processing " << processNEvents << " ouf of "
			<< nentries << " events. " << std::endl;

	tree::Photon *thisphoton = new tree::Photon();
	tree::Jet *thisjet = new tree::Jet();

	tree->Branch("photon", &photon);
	tree->Branch("jet", &jet);
	tree->Branch("met", &met, "met/F");
	tree->Branch("nVertex", &nVertex, "nVertex/I");
	tree->Branch("weigth", &weight, "weight/F");


	for (long jentry=0; jentry < processNEvents; jentry++) {
		if (jentry%reportEvery==0)
			std::cout << jentry << " / " << processNEvents << std :: endl;
		inputTree->GetEntry(jentry);

		photon.clear();
		jet.clear();

		// photons
		std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
		for(std::vector<susy::Photon>::iterator it = phoMap->second.begin();
				it != phoMap->second.end() && phoMap != event->photons.end(); it++ ) {
			if( ! it->isEB() )
				continue;
			thisphoton->pt = it->momentum.Et();
			if( thisphoton->pt < 80 )
				continue;
			thisphoton->eta = it->momentum.Eta();
			thisphoton->chargedIso = it->chargedHadronIso;
			thisphoton->neutralIso = it->neutralHadronIso;
			thisphoton->photonIso = it->photonIso;
			if ( it->r9 > 1 ) // if == 1 ?
				continue;
			thisphoton->r9 = it->r9;
			thisphoton->sigmaIetaIeta = it->sigmaIetaIeta;
			thisphoton->hadTowOverEm = it->hadTowOverEm;
			thisphoton->pixelseed = it->nPixelSeeds;
			photon.push_back( *thisphoton );
		}
		if( photon.size() == 0 )
			continue;
		std::sort( photon.begin(), photon.end(), tree::EtGreater);


		// jets
		std::map<TString,susy::PFJetCollection>::iterator pfJets_it = event->pfJets.find("ak5");
		if(pfJets_it == event->pfJets.end()){
			if(event->pfJets.size() > 0) std::cout << "JetCollection is not available!!!" << std::endl;
		} else {

			susy::PFJetCollection& jetColl = pfJets_it->second;

			for(std::vector<susy::PFJet>::iterator it = jetColl.begin();
					it != jetColl.end(); it++) {
				std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
				if (s_it == it->jecScaleFactors.end()) {
					std::cout << "JEC is not available for this jet!!!" << std::endl;
					continue;
				}
				float scale = s_it->second;
				TLorentzVector corrP4 = scale * it->momentum;

				if(std::abs(corrP4.Eta()) > 3.0) continue;
				thisjet->pt = corrP4.Et();
				thisjet->eta = corrP4.Eta();
				jet.push_back( *thisjet );
			}// for jet
		}// if, else
		if( jet.size() == 0 )
			std::cout << "error, no jets found " << std::endl;
		else
			std::sort( jet.begin(), jet.end(), tree::EtGreater);

		// met
		std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfMet");
		susy::MET* metobj = &(met_it->second);
		met = metobj->met();

		// vertices

		nVertex = event->vertices.size();
		weight = 1;

		tree->Fill();
	} // for jentry



	tree->Write();
	outFile->cd();
	outFile->Write();
	outFile->Close();
}
/*
int main(int argc, char** argv) {
	TreeWriter *tw = new TreeWriter("qcd-1000-nTuple-test.root", "myQCDTree.root");
	tw->SetProcessNEvents(10);
	tw->Loop();
}
*/
