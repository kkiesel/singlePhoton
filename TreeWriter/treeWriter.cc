#include "treeWriter.h"

using namespace std;

TreeWriter::TreeWriter( std::string inputName, std::string outputName, int loggingVerbosity_ ) {
	// read the input file
	inputTree = new TChain("susyTree");
	if (loggingVerbosity_ > 0)
		std::cout << "Add files to chain" << std::endl;
	inputTree->Add( inputName.c_str() );
	Init( outputName, loggingVerbosity_ );
}

TreeWriter::TreeWriter( TChain* inputTree_, std::string outputName, int loggingVerbosity_ ) {
	inputTree = inputTree_;
	Init( outputName, loggingVerbosity_ );
}

void TreeWriter::Init( std::string outputName, int loggingVerbosity_ ) {

	if (loggingVerbosity_ > 0)
		std::cout << "Set Branch Address of susy::Event" << std::endl;
	event = new susy::Event;
	event->setInput( *inputTree );

	// Here the number of proceeded events will be stored. For plotting, simply use L*sigma/eventNumber
	eventNumbers = new TH1F("eventNumbers", "Histogram containing number of generated events", 1, 0, 1);
	eventNumbers->GetXaxis()->SetBinLabel(1,"Number of generated events");

	// open the output file
	if (loggingVerbosity_>0)
		std::cout << "Open file " << outputName << " for writing." << std::endl;
	outFile = new TFile( outputName.c_str(), "recreate" );
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
	delete event;
	delete inputTree;
	delete outFile;
	delete tree;
}

float deltaPhi( float phi1, float phi2) {
	float result = phi1 - phi2;
	while (result > M_PI) result -= 2*M_PI;
	while (result <= -M_PI) result += 2*M_PI;
	return result;
}

// useful functions
float deltaR( const TLorentzVector& v1, const TLorentzVector& v2 ) {
	// deltaR  = sqrt ( deltaEta^2 + deltaPhi^2 )
	return sqrt(pow(v1.Eta() - v2.Eta(), 2) + pow(deltaPhi(v1.Phi(),v2.Phi()), 2) );
}

float effectiveAreaElectron( float eta ) {
	// needed by calculating the isolation for electrons
	// see https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection
	// only for Delta R = 0.3 on 2012 Data
	eta = fabs( eta );
	float ea;
	if( eta < 1.0 ) ea = 0.13;
	else if( eta < 1.479 ) ea = 0.14;
	else if( eta < 2.0 ) ea = 0.07;
	else if( eta < 2.2 ) ea = 0.09;
	else if( eta < 2.3 ) ea = 0.11;
	else if( eta < 2.4 ) ea = 0.11;
	else ea = 0.14;
	return ea;
}

// correct iso, see https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
float chargedHadronIso_corrected(const susy::Photon& gamma, float rho) {
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

float neutralHadronIso_corrected(const susy::Photon& gamma, float rho) {
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

float photonIso_corrected(const susy::Photon& gamma, float rho) {
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

float d0correction( const susy::Electron& electron, const susy::Event& event ) {
	// copied from Brian Francis
	TVector3 beamspot = event.vertices[0].position;
	susy::Track track = event.tracks[electron.gsfTrackIndex];
	float d0 = track.d0() - beamspot.X()*sin(track.phi()) + beamspot.Y()*cos(track.phi());
	return d0;
}

float dZcorrection( const susy::Electron& electron, const susy::Event& event ) {
	// copied from Brian Francis
	TVector3 beamspot = event.vertices[0].position;
	susy::Track track = event.tracks[electron.gsfTrackIndex];

	if(track.momentum.Pt() == 0.) return 1.e6;
	float dz = (track.vertex.Z() - beamspot.Z()) - ((track.vertex.X() - beamspot.X())*track.momentum.Px() + (track.vertex.Y() - beamspot.Y())*track.momentum.Py()) / track.momentum.Pt() * (track.momentum.Pz() / track.momentum.Pt());
	return dz;
}

float getPtFromMatchedJet( const susy::Photon& myPhoton, const susy::PFJetCollection& jetColl, int loggingVerbosity = 0 ) {
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

	for(std::vector<susy::PFJet>::const_iterator it = jetColl.begin();
			it != jetColl.end(); ++it) {
		float scale = 1.;
		std::map<TString,Float_t>::const_iterator s_it = it->jecScaleFactors.find("L2L3");
		if (s_it == it->jecScaleFactors.end()) {
			std::cout << "JEC is not available for this jet!!!" << std::endl;
			continue;
		} else {
			scale = s_it->second;
		}
		TLorentzVector corrP4 = scale * it->momentum;
		float deltaR_ = deltaR(myPhoton.momentum, corrP4 );
		if (deltaR_ > 0.3) continue;
		if( loggingVerbosity > 2 )
			std::cout << " pT_jet / pT_gamma = " << it->momentum.Et() / myPhoton.momentum.Et() << std::endl;
		nearJets.push_back( *it );
	}// for jet

	if ( nearJets.size() == 0 ) {
		if( loggingVerbosity > 1 )
			std::cout << "No jet with deltaR < .3 found, do not change photon_pt" << std::endl;
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
	/**
	 * \brief Loops over input chain and fills tree
	 *
	 * This is the major function of treeWriter, which initialize the output, loops
	 * over all events and fill the tree. In the end, the tree is saved to the
	 * output File
	 */

	// here the event loop is implemented and the tree is filled
	if (inputTree == 0) return;

	// get number of events to be proceeded
	Long64_t nentries = inputTree->GetEntries();
	// store them in histo
	eventNumbers->Fill( "Number of generated events", nentries );
	if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

	if( loggingVerbosity > 0 )
		std::cout << "Processing " << processNEvents << " ouf of "
			<< nentries << " events. " << std::endl;

	tree->Branch("photon", &photon);
	tree->Branch("jet", &jet);
	tree->Branch("electron", &electron);
	tree->Branch("muon", &muon);
	tree->Branch("met", &met, "met/F");
	tree->Branch("metPhi", &met_phi, "metPhi/F");
	tree->Branch("type1met", &type1met, "type1met/F");
	tree->Branch("type1metPhi", &type1met_phi, "type1metPhi/F");
	tree->Branch("ht", &ht, "ht/F");
	tree->Branch("nVertex", &nVertex, "nVertex/I");
	tree->Branch("weight", &weight, "weight/D");
	tree->Branch("genElectron", &genElectron);
	tree->Branch("genPhoton", &genPhoton);

	for (long jentry=0; jentry < processNEvents; ++jentry) {
		if ( loggingVerbosity>1 || jentry%reportEvery==0 )
			std::cout << jentry << " / " << processNEvents << std :: endl;
		inputTree->LoadTree( jentry );
		event->getEntry(jentry);

		photon.clear();
		jet.clear();
		electron.clear();
		muon.clear();
		genElectron.clear();
		genPhoton.clear();
		ht = 0;

		// weights
		if (pileupHisto == 0) {
			weight = 1.;
		} else {
			float trueNumInteractions = -1;
			for( susy::PUSummaryInfoCollection::const_iterator iBX = event->pu.begin();
					iBX != event->pu.end() && trueNumInteractions < 0; ++iBX) {
				if (iBX->BX == 0)
					trueNumInteractions = iBX->trueNumInteractions;
			}
			weight = pileupHisto->GetBinContent( pileupHisto->FindBin( trueNumInteractions ) );
		}

		// get ak5 jets
		std::vector<susy::PFJet> jetVector = event->pfJets["ak5"];

		// photons
		std::vector<susy::Photon> photonVector = event->photons["photons"];

		for(std::vector<susy::Photon>::iterator it = photonVector.begin();
				it != photonVector.end(); ++it ) {
			if( !(it->isEE() || it->isEB()) && it->momentum.Pt()<20 && it->isEBEtaGap() && it->isEBPhiGap() && it->isEERingGap() && it->isEEDeeGap() && it->isEBEEGap() && skim )
				continue;
			tree::Photon thisphoton;

			thisphoton.chargedIso = chargedHadronIso_corrected(*it, event->rho25);
			thisphoton.neutralIso = neutralHadronIso_corrected(*it, event->rho25);
			thisphoton.photonIso = photonIso_corrected(*it, event->rho25);

			bool loose_photon_barrel = it->isEB()
				&& it->hadTowOverEm<0.05
				&& it->sigmaIetaIeta<0.012
				&& thisphoton.chargedIso<2.6
				&& thisphoton.neutralIso<3.5+0.04*thisphoton.pt
				&& thisphoton.photonIso<1.3+0.005*thisphoton.pt;

			bool loose_photon_endcap = it->isEE()
				&& it->hadTowOverEm<0.05
				&& it->sigmaIetaIeta<0.034
				&& thisphoton.chargedIso<2.3
				&& thisphoton.neutralIso<2.9+0.04*thisphoton.pt;

			thisphoton.ptJet = getPtFromMatchedJet( *it, jetVector, loggingVerbosity );
			if(!(loose_photon_endcap || loose_photon_barrel || thisphoton.ptJet > 75 ) && skim )
				continue;
			thisphoton.pt = it->momentum.Pt();
			thisphoton.eta = it->momentum.Eta();
			thisphoton.phi = it->momentum.Phi();
			thisphoton.r9 = it->r9;
			thisphoton.sigmaIetaIeta = it->sigmaIetaIeta;
			thisphoton.hadTowOverEm = it->hadTowOverEm;
			thisphoton.pixelseed = it->nPixelSeeds;
			thisphoton.conversionSafeVeto = it->passelectronveto;
			thisphoton.genInformation = 0;
			photon.push_back( thisphoton );
			if( loggingVerbosity > 2 )
				std::cout << " p_T, gamma = " << thisphoton.pt << std::endl;
		}

		if( photon.size() == 0 && skim )
			continue;
		std::sort( photon.begin(), photon.end(), tree::EtGreater);
		if( loggingVerbosity > 1 )
			std::cout << "Found " << photon.size() << " photons" << std::endl;

		// jets


		for(std::vector<susy::PFJet>::iterator it = jetVector.begin();
				it != jetVector.end(); ++it) {
			tree::Jet thisjet;

			// scale with JEC
			float scale = 1.;
			if(it->jecScaleFactors.count("L2L3") == 0)
				std::cout << "ERROR: JEC is not available for this jet" << std::endl;
			else
				scale = it->jecScaleFactors.find("L2L3")->second;
			TLorentzVector corrP4 = scale * it->momentum;

			// Calculate HT.
			// The definiton differs from the saved jet, since trigger is described better
			if( std::abs( corrP4.Eta() ) < 3 && corrP4.Pt() > 40 )
				ht += thisjet.pt;

			if(std::abs(corrP4.Eta()) > 2.6 && skim ) continue;
			if(corrP4.Pt() < 30 && skim ) continue;
			thisjet.pt = corrP4.Pt();
			thisjet.eta = corrP4.Eta();
			thisjet.phi = corrP4.Phi();
			thisjet.bCSV = it->bTagDiscriminators[susy::kCSV];
			// jet composition
			thisjet.chargedHadronEnergy = it->chargedHadronEnergy;
			thisjet.neutralHadronEnergy = it->neutralHadronEnergy;
			thisjet.photonEnergy = it->photonEnergy;
			thisjet.electronEnergy = it->electronEnergy;
			thisjet.muonEnergy = it->muonEnergy;
			thisjet.HFHadronEnergy = it->HFHadronEnergy;
			thisjet.HFEMEnergy = it->HFEMEnergy;
			thisjet.chargedEmEnergy = it->chargedEmEnergy;
			thisjet.chargedMuEnergy = it->chargedMuEnergy;
			thisjet.neutralEmEnergy = it->neutralEmEnergy;

			if( loggingVerbosity > 2 )
				std::cout << " p_T, jet = " << thisjet.pt << std::endl;

			jet.push_back( thisjet );
		}// for jet
		if( jet.size() < 2 && skim )
			continue;
		std::sort( jet.begin(), jet.end(), tree::EtGreater);
		if( loggingVerbosity > 1 )
			std::cout << "Found " << jet.size() << " jets" << std::endl;

		if( ht < 450 && skim)
			continue;



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
		std::vector<susy::Electron> eVector = event->electrons["gsfElectrons"];
		for(std::vector<susy::Electron>::iterator it = eVector.begin(); it < eVector.end(); ++it) {
			tree::Particle thiselectron;
			if( loggingVerbosity > 2 )
				cout << " electron pt = " << it->momentum.Pt() << endl;
			// for cuts see https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
			// use veto electrons
			if( it->momentum.Pt() < 20  || it->momentum.Pt() > 1e6 )
				continue; // spike rejection
			float iso = ( it->chargedHadronIso + max(it->neutralHadronIso+it->photonIso - effectiveAreaElectron(it->momentum.Eta())*event->rho25, (float)0. )
						) / it->momentum.Pt();
			float d0 = d0correction( *it, *event );
			float dZ = std::abs( dZcorrection( *it, *event ) );
			if ( it->isEB() ){
				if ( fabs(it->deltaEtaSuperClusterTrackAtVtx) > 0.007
						|| fabs(it->deltaPhiSuperClusterTrackAtVtx) > 0.8
						|| it->sigmaIetaIeta > 0.01
						|| it->hcalOverEcalBc > 0.15
						|| d0 > 0.04
						|| dZ > 0.2
						|| iso > 0.15 )
					continue;
				}
			else if( it->isEE() ) {
				if ( fabs(it->deltaEtaSuperClusterTrackAtVtx) > 0.01
						|| fabs(it->deltaPhiSuperClusterTrackAtVtx) > 0.7
						|| it->sigmaIetaIeta > 0.03
						|| d0 > 0.04
						|| dZ > 0.2
						|| iso > 0.15 )
					continue;
				}
			else // not in barrel nor in endcap
				continue;

			thiselectron.pt = it->momentum.Pt();
			if( loggingVerbosity > 2 )
				std::cout << " p_T, electron = " << it->momentum.Et() << std::endl;
			thiselectron.eta = it->momentum.Eta();
			thiselectron.phi = it->momentum.Phi();
			electron.push_back( thiselectron );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << electron.size() << " electrons" << std::endl;

		// muons
		tree::Particle thismuon;
		std::vector<susy::Muon> mVector = event->muons["muons"];
		for( std::vector<susy::Muon>::iterator it = mVector.begin(); it != mVector.end(); ++it) {
			if( !( it->isPFMuon() && ( it->isGlobalMuon() || it->isTrackerMuon() ) ) )
				continue; // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
			thismuon.pt = it->momentum.Et();
			if( thismuon.pt < 20 )
				continue;
			thismuon.eta = it->momentum.Eta();
			thismuon.phi = it->momentum.Phi();
			muon.push_back( thismuon );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << muon.size() << " muons" << std::endl;

		// vertices
		nVertex = event->vertices.size();

		tree::Particle thisGenParticle;
		for( std::vector<susy::Particle>::iterator it = event->genParticles.begin(); it != event->genParticles.end(); ++it ) {
			if( it->momentum.Pt() < 20 ) continue;
			thisGenParticle.pt = it->momentum.Pt();
			thisGenParticle.eta = it->momentum.Eta();
			thisGenParticle.phi = it->momentum.Phi();
			switch( std::abs(it->pdgId) ) {
				case 22: // photon
					genPhoton.push_back( thisGenParticle );
					break;
				case 11: // electron
					// Demand a W boson as mother particle of electron
					if( abs(event->genParticles[it->motherIndex].pdgId) == 24 )
						genElectron.push_back( thisGenParticle );
					break;
			}
		}

		tree->Fill();
	} // for jentry


	outFile->cd();
	eventNumbers->Write();
	tree->Write();
	outFile->Write();
	outFile->Close();
}

