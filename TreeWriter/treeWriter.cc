#include<iostream>
#include<math.h>

#include "TSystem.h"

#include "treeWriter.h"

using namespace std;

TreeWriter::TreeWriter(TString inputName, TString outputName, int loggingVerbosity_ ) {
	// read the input file
	inputTree = new TChain("susyTree");
	if (loggingVerbosity_ > 0)
		std::cout << "Add files to chain" << std::endl;
	inputTree->Add( inputName );

	if (loggingVerbosity_ > 0)
		std::cout << "Set Branch Address of susy::Event" << std::endl;
	event = new susy::Event;
	inputTree->SetBranchAddress("susyEvent", &event);

	// open the output file
	outFile = new TFile( outputName, "recreate" );
	tree = new TTree("susyTree","Tree for single photon analysis");

	// set default parameter
	processNEvents = -1;
	reportEvery = 1000;
	loggingVerbosity = loggingVerbosity_;
	skim = true;

}

TreeWriter::~TreeWriter() {
	if (!inputTree) return;
	delete inputTree->GetCurrentFile();
}

float TreeWriter::deltaR( TLorentzVector v1, TLorentzVector v2 ) {
	return sqrt(pow(v1.Eta() - v2.Eta(), 2) + pow(v1.Phi() - v2.Phi(), 2) );
}

float TreeWriter::getPtFromMatchedJet( susy::Photon myPhoton, susy::Event myEvent ) {
	/**
	 * \brief Takes jet p_T as photon p_T
	 *
	 * At first all jets with DeltaR < 0.3 (isolation cone) are searched.
	 * If several jets are found, take the one with the minimal pt difference
	 * compared to the photon. If no such jets are found, keep the photon_pt
	 * TODO: remove photon matched jet from jet-selection?
	 */
	std::vector<susy::PFJet> nearJets;
	nearJets.clear();

	std::map<TString,susy::PFJetCollection>::iterator pfJets_it = myEvent.pfJets.find("ak5");
	if(pfJets_it == myEvent.pfJets.end()){
		if(myEvent.pfJets.size() > 0) std::cout << "JetCollection is not available!!!" << std::endl;
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
			float deltaR_ = deltaR(myPhoton.momentum, corrP4 );
			if (deltaR_ > 0.3) continue;
			if( loggingVerbosity > 0 )
				std::cout << "gamma pt jet matching factor = " << it->momentum.Et() / myPhoton.momentum.Et() << std::endl;
			nearJets.push_back( *it );
		}// for jet
	}// if, else

	if ( nearJets.size() == 0 ) {
		//std::cout << "No jet with deltaR < .3 found, do not change photon_pt" << std::endl;
		return myPhoton.momentum.Et();
	}

	float pt = 0;
	float minPtDifferenz = 1E20; // should be very high
	for( std::vector<susy::PFJet>::iterator it = nearJets.begin(), jetEnd = nearJets.end();
			it != jetEnd; it++ ) {
		float ptDiff = fabs(myPhoton.momentum.Et() - it->momentum.Et());
		if (  ptDiff < minPtDifferenz ) {
			minPtDifferenz = ptDiff;
			pt = it->momentum.Et();
		}
	}

	// testing
	if( nearJets.size() > 1 && loggingVerbosity > 0 )
		std::cout << "There are several jets matching to this photon. "
					<< "Please check if jet-matching is correct." << std::endl;
	return pt;
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
	tree->Branch("ht", &ht, "ht/F");
	tree->Branch("nVertex", &nVertex, "nVertex/I");
	tree->Branch("nElectron", &nElectron, "nElectron/I");


	for (long jentry=0; jentry < processNEvents; jentry++) {
		if ( loggingVerbosity>0 && jentry%reportEvery==0 )
			std::cout << jentry << " / " << processNEvents << std :: endl;
		inputTree->GetEntry(jentry);

		photon.clear();
		jet.clear();
		ht = 0;

		// photons
		if( loggingVerbosity > 1 )
			std::cout << "Process photons" << std::endl;
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
			ht += thisphoton->pt;
			if( loggingVerbosity > 2 )
				std::cout << " p_T, gamma = " << thisphoton->pt << std::endl;
		}

		if( photon.size() == 0 && skim )
			continue;
		std::sort( photon.begin(), photon.end(), tree::EtGreater);
		if( loggingVerbosity > 1 )
			std::cout << "Found " << photon.size() << " photons" << std::endl;

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
				ht += thisjet->pt;
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
		//std::vector<susy::Electron> eVector = event->electrons["gsfElectrons"];
		//nElectron = eVector.size();

		// vertices

		nVertex = event->vertices.size();

		if( ht < 450 && skim)
			continue;

		tree->Fill();
	} // for jentry



	tree->Write();
	outFile->cd();
	outFile->Write();
	outFile->Close();
}

