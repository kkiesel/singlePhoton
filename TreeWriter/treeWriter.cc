#include<iostream>
#include<math.h>
#include<string>

#include "TSystem.h"

#include "treeWriter.h"

using namespace std;

TreeWriter::TreeWriter( TString inputName, TString outputName, int loggingVerbosity_ ) {
	// read the input file
	inputTree = new TChain("susyTree");
	if (loggingVerbosity_ > 0)
		std::cout << "Add files to chain" << std::endl;
	inputTree->Add( inputName );
	Init( outputName, loggingVerbosity_ );
}

TreeWriter::TreeWriter( TChain* inputTree_, TString outputName, int loggingVerbosity_ ) {
	inputTree = inputTree_;
	Init( outputName, loggingVerbosity_ );
}


void TreeWriter::Init( TString outputName, int loggingVerbosity_ ) {

	if (loggingVerbosity_ > 0)
		std::cout << "Set Branch Address of susy::Event" << std::endl;
	event = new susy::Event;
	inputTree->SetBranchAddress("susyEvent", &event);

	// open the output file
	if (loggingVerbosity_>0)
		std::cout << "Open file " << outputName << " for writing." << std::endl;
	outFile = new TFile( outputName, "recreate" );
	tree = new TTree("susyTree","Tree for single photon analysis");

	// set default parameter
	processNEvents = -1;
	reportEvery = 1000;
	loggingVerbosity = loggingVerbosity_;
	skim = true;
	pileupHisto = 0;
}

void TreeWriter::PileUpWeightFile( string pileupFileName ) {
	TFile *puFile = new TFile( pileupFileName.c_str() );
	pileupHisto = (TH1F*) puFile->Get("pileup");
}

TreeWriter::~TreeWriter() {
	if (pileupHisto != 0 )
		delete pileupHisto;
	inputTree->GetCurrentFile()->Close();
}

// useful functions
float TreeWriter::deltaR( TLorentzVector v1, TLorentzVector v2 ) {
	return sqrt(pow(v1.Eta() - v2.Eta(), 2) + pow(v1.Phi() - v2.Phi(), 2) );
}

// correct iso, see https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
float chargedHadronIso_corrected(susy::Photon gamma, float rho) {
  float eta = fabs(gamma.caloPosition.Eta());
  float ea;

  if(eta < 1.0) ea = 0.012;
  else if(eta < 1.479) ea = 0.010;
  else if(eta < 2.0) ea = 0.014;
  else if(eta < 2.2) ea = 0.012;
  else if(eta < 2.3) ea = 0.016;
  else if(eta < 2.4) ea = 0.020;
  else ea = 0.012;

  float iso = gamma.chargedHadronIso;
  iso = max(iso - rho*ea, (float)0.);

  return iso;
}

float neutralHadronIso_corrected(susy::Photon gamma, float rho) {
  float eta = fabs(gamma.caloPosition.Eta());
  float ea;

  if(eta < 1.0) ea = 0.030;
  else if(eta < 1.479) ea = 0.057;
  else if(eta < 2.0) ea = 0.039;
  else if(eta < 2.2) ea = 0.015;
  else if(eta < 2.3) ea = 0.024;
  else if(eta < 2.4) ea = 0.039;
  else ea = 0.072;

  float iso = gamma.neutralHadronIso;
  iso = max(iso - rho*ea, (float)0.);

  return iso;
}

float photonIso_corrected(susy::Photon gamma, float rho) {
  float eta = fabs(gamma.caloPosition.Eta());
  float ea;

  if(eta < 1.0) ea = 0.148;
  else if(eta < 1.479) ea = 0.130;
  else if(eta < 2.0) ea = 0.112;
  else if(eta < 2.2) ea = 0.216;
  else if(eta < 2.3) ea = 0.262;
  else if(eta < 2.4) ea = 0.260;
  else ea = 0.266;

  float iso = gamma.photonIso;
  iso = max(iso - rho*ea, (float)0.);

  return iso;
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
				it != jetColl.end(); ++it) {
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
			it != jetEnd; ++it ) {
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
	//tree->Branch("electron", &electron);
	tree->Branch("muon", &muon);
	tree->Branch("met", &met, "met/F");
	tree->Branch("metPhi", &met_phi, "metPhi/F");
	tree->Branch("type1met", &type1met, "type1met/F");
	tree->Branch("type1metPhi", &type1met_phi, "type1metPhi/F");
	tree->Branch("ht", &ht, "ht/F");
	tree->Branch("nVertex", &nVertex, "nVertex/I");
	tree->Branch("pu_weight", &pu_weight, "pu_weight/F");


	for (long jentry=0; jentry < processNEvents; ++jentry) {
		if ( loggingVerbosity>0 && jentry%reportEvery==0 )
			std::cout << jentry << " / " << processNEvents << std :: endl;
		inputTree->GetEntry(jentry);

		photon.clear();
		jet.clear();
		electron.clear();
		muon.clear();
		ht = 0;

		// weights
		if (pileupHisto == 0) {
			pu_weight = 1.;
		} else {
			float trueNumInteractions = -1;
			for( susy::PUSummaryInfoCollection::const_iterator iBX = event->pu.begin();
					iBX != event->pu.end() && trueNumInteractions < 0; ++iBX) {
				if (iBX->BX == 0)
					trueNumInteractions = iBX->trueNumInteractions;
			}
			pu_weight = pileupHisto->GetBinContent( pileupHisto->FindBin( trueNumInteractions ) );
		}


		// photons
		if( loggingVerbosity > 1 )
			std::cout << "Process photons" << std::endl;
		std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
		for(std::vector<susy::Photon>::iterator it = phoMap->second.begin();
				it != phoMap->second.end() && phoMap != event->photons.end(); ++it ) {
			if( ! it->isEB() && skim )
				continue;
			thisphoton->pt = getPtFromMatchedJet( *it, *event );
			if( thisphoton->pt < 80 && skim )
				continue;
			thisphoton->eta = it->momentum.Eta();
			thisphoton->phi = it->momentum.Phi();
			thisphoton->chargedIso = chargedHadronIso_corrected(*it, event->rho25);
			thisphoton->neutralIso = neutralHadronIso_corrected(*it, event->rho25);
			thisphoton->photonIso = photonIso_corrected(*it, event->rho25);
			thisphoton->r9 = it->r9;
			thisphoton->sigmaIetaIeta = it->sigmaIetaIeta;
			thisphoton->hadTowOverEm = it->hadTowOverEm;
			thisphoton->pixelseed = it->nPixelSeeds;
			thisphoton->conversionSafeVeto = it->passelectronveto;
			photon.push_back( *thisphoton );
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
					it != jetColl.end(); ++it) {
				std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
				if (s_it == it->jecScaleFactors.end()) {
					std::cout << "JEC is not available for this jet!!!" << std::endl;
					continue;
				}
				float scale = s_it->second;
				TLorentzVector corrP4 = scale * it->momentum;

				if(std::abs(corrP4.Eta()) > 3.0 && skim ) continue;
				if(corrP4.Et() < 30 && skim ) continue;
				thisjet->pt = corrP4.Et();
				thisjet->eta = corrP4.Eta();
				thisjet->phi = corrP4.Phi();
				thisjet->bCSV = it->bTagDiscriminators[susy::kCSV];
				if( loggingVerbosity > 2 )
					std::cout << " p_T, jet = " << thisjet->pt << std::endl;

				jet.push_back( *thisjet );
				ht += thisjet->pt;
			}// for jet
		}// if, else
		if( jet.size() < 2 && skim )
			continue;
		std::sort( jet.begin(), jet.end(), tree::EtGreater);
		if( loggingVerbosity > 1 )
			std::cout << "Found " << jet.size() << " jets" << std::endl;


		// met
		std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfMet");
		susy::MET* metobj = &(met_it->second);
		met = metobj->met();
		met_phi = metobj->mEt.Phi();
		if( loggingVerbosity > 2 )
			std::cout << " met = " << met << std::endl;

		std::map<TString, susy::MET>::iterator type1met_it = event->metMap.find("pfType1CorrectedMet");
		susy::MET* type1metobj = &(type1met_it->second);
		type1met = type1metobj->met();
		type1met_phi = type1metobj->mEt.Phi();
		if( loggingVerbosity > 2 )
			std::cout << " type1met = " << type1met << std::endl;

		// electrons
		/*
		tree::Particle* thiselectron = new tree::Particle();
		map<TString, vector<susy::Electron> >::iterator eleMap = event->electrons.find("gsfElectrons");
		if(eleMap == event->electrons.end() && loggingVerbosity > 0) {
			cout << "gsfElectrons not found!" << endl;
		} else {
			for(vector<susy::Electron>::iterator it = eleMap->second.begin(); it < eleMap->second.end(); ++it) {
				thiselectron->pt = it->momentum.Et();
				if( thiselectron->pt < 20 )
					continue;
				if( loggingVerbosity > 2 )
					std::cout << " p_T, electron = " << it->momentum.Et() << std::endl;
				thiselectron->eta = it->momentum.Eta();
				thiselectron->phi = it->momentum.Phi();
				electron.push_back( *thiselectron );
			}
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << electron.size() << " electrons" << std::endl;
		*/

		// this seems not to work yet, where is the bug?
		/*
		std::vector<susy::Electron> eVector = event->electrons["gsfElectronsx"];
		for( std::vector<susy::Electron>::iterator it = eVector.begin(); it != eVector.end(); ++it) {
			thiselectron->pt = it->momentum.Et();
			if( thiselectron->pt < 20 )
				continue;
			if( loggingVerbosity > 2 )
				std::cout << " p_T, electron = " << it->momentum.Et() << std::endl;
			thiselectron->eta = it->momentum.Eta();
			thiselectron->phi = it->momentum.Phi();
			electron.push_back( *thiselectron );
		}
		*/

		// muons
		std::vector<susy::Muon> mVector = event->muons["muons"];
		tree::Particle* thismuon = new tree::Particle();
		for( std::vector<susy::Muon>::iterator it = mVector.begin(); it != mVector.end(); ++it) {
			if( !( it->isPFMuon() && ( it->isGlobalMuon() || it->isTrackerMuon() ) ) )
				continue; // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
			thismuon->pt = it->momentum.Et();
			if( thismuon->pt < 20 )
				continue;
			thismuon->eta = it->momentum.Eta();
			thismuon->phi = it->momentum.Phi();
			muon.push_back( *thismuon );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << muon.size() << " muons" << std::endl;


		// vertices
		nVertex = event->vertices.size();

		if( ht < 450 && skim)
			continue;

		tree->Fill();
	} // for jentry


	outFile->cd();
	tree->Write();
	outFile->Write();
	outFile->Close();
}

