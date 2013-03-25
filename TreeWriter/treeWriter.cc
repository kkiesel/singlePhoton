#include<iostream>
#include<math.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "SusyEvent.h"
#include "TreeObjects.h"


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

		float getPtFromMatchedJet( susy::Photon, susy::Event );
		float deltaR( TLorentzVector, TLorentzVector );

	private:
		int processNEvents; // number of events to be processed
		int reportEvery;
		int loggingVerbosity;

		// variables which will be stored in the tree
		std::vector<tree::Photon> photon;
		std::vector<tree::Jet> jet;
		float met;
		int nVertex;
		int nElectron;
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

float TreeWriter::deltaR( TLorentzVector v1, TLorentzVector v2 ) {
	return sqrt(pow(v1.Eta() - v2.Eta(), 2) + pow(v1.Phi() - v2.Phi(), 2) );
}

float TreeWriter::getPtFromMatchedJet( susy::Photon photon, susy::Event event ) {
	/**
	 * \brief Takes jet p_T as photon p_T
	 *
	 * At first all jets with DeltaR < 0.3 (isolation cone) are searched.
	 * If several jets are found, take the one with the minimal pt difference
	 * compared to the photon. If no such jets are found, return 0
	 * TODO: remove photon matched jet from jet-selection?
	 */
	std::vector<susy::PFJet> nearJets;
	nearJets.clear();

	std::map<TString,susy::PFJetCollection>::iterator pfJets_it = event.pfJets.find("ak5");
	if(pfJets_it == event.pfJets.end()){
		if(event.pfJets.size() > 0) std::cout << "JetCollection is not available!!!" << std::endl;
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
			float deltaR_ = deltaR(photon.momentum, corrP4 );
			if (deltaR_ > 0.3) continue;
			nearJets.push_back( *it );
		}// for jet
	}// if, else

	if ( nearJets.size() == 0 ) {
		std::cout << "No jet with ΔR < .3 found, set p_T to 0" << std::endl;
		return 0;
	}

	float pt = 0;
	float minPtDifferenz = 1E20; // should be very high
	for( std::vector<susy::PFJet>::iterator it = nearJets.begin(), jetEnd = nearJets.end();
			it != jetEnd; it++ ) {
		float ptDiff = fabs(photon.momentum.Et() - it->momentum.Et());
		if (  ptDiff < minPtDifferenz ) {
			minPtDifferenz = ptDiff;
			pt = it->momentum.Et();
		}
	}

	// testing
	if( nearJets.size() > 1 )
		std::cout << "There are several jets matching to this photon. "
					<< "Please check if the upper assumtion is reasonable" << std::endl;
	return pt;
}


void TreeWriter::Loop() {
	// only for testing
	bool skim = true;

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
	tree->Branch("nElectron", &nElectron, "nElectron/I");


	for (long jentry=0; jentry < processNEvents; jentry++) {
		if ( loggingVerbosity>0 && jentry%reportEvery==0 )
			std::cout << jentry << " / " << processNEvents << std :: endl;
		inputTree->GetEntry(jentry);

		photon.clear();
		jet.clear();

		// photons
		std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
		for(std::vector<susy::Photon>::iterator it = phoMap->second.begin();
				it != phoMap->second.end() && phoMap != event->photons.end(); it++ ) {
			if( ! it->isEB() && skim )
				continue;
			thisphoton->pt = getPtFromMatchedJet( *it, *event );
			if( thisphoton->pt < 80 && skim )
				continue;
			thisphoton->eta = it->momentum.Eta();
			thisphoton->phi = it->momentum.Phi();
			thisphoton->chargedIso = it->chargedHadronIso;
			thisphoton->neutralIso = it->neutralHadronIso;
			thisphoton->photonIso = it->photonIso;
			if ( it->r9 > 1 && skim ) // if == 1 ?
				continue;
			thisphoton->r9 = it->r9;
			thisphoton->sigmaIetaIeta = it->sigmaIetaIeta;
			thisphoton->hadTowOverEm = it->hadTowOverEm;
			thisphoton->pixelseed = it->nPixelSeeds;
			photon.push_back( *thisphoton );
		}
		if( photon.size() == 0 && skim )
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

				if(std::abs(corrP4.Eta()) > 3.0 && skim ) continue;
				thisjet->pt = corrP4.Et();
				thisjet->eta = corrP4.Eta();
				thisjet->phi = corrP4.Phi();
				jet.push_back( *thisjet );
			}// for jet
		}// if, else
		if( jet.size() < 2 && skim )
			continue;
		std::sort( jet.begin(), jet.end(), tree::EtGreater);

		// met
		std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfMet");
		susy::MET* metobj = &(met_it->second);
		met = metobj->met();

		// electrons
		std::vector<susy::Electron> eVector = event->electrons["gsfElectrons"];
		nElectron = eVector.size();

		// vertices

		nVertex = event->vertices.size();
		weight = 1;

		// bjets

		tree->Fill();
	} // for jentry



	tree->Write();
	outFile->cd();
	outFile->Write();
	outFile->Close();
}

int main(int argc, char** argv) {
	TreeWriter *tw = new TreeWriter("../root-files/susyEvents_qcd_1000-inf_part.root", "myTree.root");
	tw->SetProcessNEvents(-1);
	tw->Loop();
}

