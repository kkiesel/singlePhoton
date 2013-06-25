#include "treeWriter.h"

float deltaPhi( float phi1, float phi2) {
	/** Delta Phi is computed for zylindical coordinates in the CMS system.
	 */
	float result = phi1 - phi2;
	while (result > M_PI) result -= 2*M_PI;
	while (result <= -M_PI) result += 2*M_PI;
	return result;
}

float deltaR( const TLorentzVector& v1, const TLorentzVector& v2 ) {
	// deltaR  = sqrt ( deltaEta^2 + deltaPhi^2 )
	return sqrt(pow(v1.Eta() - v2.Eta(), 2) + pow(deltaPhi(v1.Phi(),v2.Phi()), 2) );
}

float deltaR( const susy::PFJet& v1, const tree::Particle& v2 ) {
	// deltaR  = sqrt ( deltaEta^2 + deltaPhi^2 )
	return sqrt(pow(v1.momentum.Eta() - v2.eta, 2) + pow(deltaPhi(v1.momentum.Phi(),v2.phi), 2) );
}

float effectiveAreaElectron( float eta ) {
	/** Returns the effective area for the isolation criteria for electrons.
	 * See https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection
	 * only for Delta R = 0.3 on 2012 Data
	 */
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

float chargedHadronIso_corrected(const susy::Photon& gamma, float rho) {
	/** Correct isolation for photons,
	 * see https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
	 */
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
	iso = std::max(iso - rho*ea, (float)0.);

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
	iso = std::max(iso - rho*ea, (float)0.);

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
	iso = std::max(iso - rho*ea, (float)0.);

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

bool isAdjacentToLepton( const susy::PFJet& jet, const std::vector<tree::Particle>& leptons, float deltaR_ = 0.3 ) {
	/** Leptons near the jet are searched.
	 *
	 * Returns true if a lepton is found in a certain radius near the jet.
	 */
	bool foundLepton = false;
	for(std::vector<tree::Particle>::const_iterator lepton = leptons.begin(); lepton != leptons.end(); ++lepton) {
		if (deltaR(jet, *lepton ) < deltaR_)
			foundLepton = true;
	}
	return foundLepton;
}

bool isVetoElectron( const susy::Electron& electron, const susy::Event& event, const int loggingVerbosity ) {
	/** Definition of veto working point for electrons.
	 *
	 * See https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
	 * for more information.
	 */
	if( electron.momentum.Pt() > 1e6 )
		return false; // spike rejection
	float iso = ( electron.chargedHadronIso +
		std::max(electron.neutralHadronIso+electron.photonIso -
		effectiveAreaElectron(electron.momentum.Eta())*event.rho25, (float)0. ))
		/ electron.momentum.Pt();
	float d0 = d0correction( electron, event );
	float dZ = std::abs( dZcorrection( electron, event ) );
	bool isElectron = false;
	float eta = std::abs(electron.momentum.Eta());
	isElectron  = (
		eta < 1.442
			&& ( fabs(electron.deltaEtaSuperClusterTrackAtVtx) < 0.007
				|| fabs(electron.deltaPhiSuperClusterTrackAtVtx) < 0.8
				|| electron.sigmaIetaIeta < 0.01
				|| electron.hcalOverEcalBc < 0.15
				|| d0 < 0.04
				|| dZ < 0.2
				|| iso < 0.15 )
		)||( (1.566 < eta && eta < 2.5)
			&& ( fabs(electron.deltaEtaSuperClusterTrackAtVtx) < 0.01
				|| fabs(electron.deltaPhiSuperClusterTrackAtVtx) < 0.7
				|| electron.sigmaIetaIeta < 0.03
				|| d0 < 0.04
				|| dZ < 0.2
				|| iso < 0.15 )
		);
	return isElectron;
}

float getPtFromMatchedJet( const susy::Photon& myPhoton, const susy::PFJetCollection& jetColl, int loggingVerbosity = 0 ) {
	/**
	 * \brief Takes jet p_T as photon p_T
	 *
	 * At first all jets with DeltaR < 0.3 (isolation cone) are searched.
	 * If several jets are found, take the one with the minimal pt difference
	 * compared to the photon. If no such jets are found, keep the photon_pt
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
	return pt;
}

///////////////////////////////////////////////////////////////////////////////
// Here the class implementation begins ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

TreeWriter::TreeWriter( std::string inputName, std::string outputName, int loggingVerbosity_ ) {
	/** Constructor function for a single input file.
	 *
	 * The tree from the file will be read.
	 * inputName: points to a root file with a "susyTree"
	 */
	inputTree = new TChain("susyTree");
	if (loggingVerbosity_ > 0)
		std::cout << "Add files to chain" << std::endl;
	inputTree->Add( inputName.c_str() );
	Init( outputName, loggingVerbosity_ );
}

TreeWriter::TreeWriter( TChain* inputTree_, std::string outputName, int loggingVerbosity_ ) {
	/** Constructor function for a TTree as input.
	 */
	inputTree = inputTree_;
	Init( outputName, loggingVerbosity_ );
}

TreeWriter::~TreeWriter() {
	/** Deconstructor
	 *
	 * Event has to be deleted before the tree.
	 */
	if (pileupHisto != 0 )
		delete pileupHisto;
	inputTree->GetCurrentFile()->Close();
	delete event;
	delete inputTree;
	delete outFile;
	delete tree;
}

void TreeWriter::Init( std::string outputName, int loggingVerbosity_ ) {
	/** Initialize the Object
	 *
	 * This function is called by each constructor function. Some basic setting
	 * are done.
	 */
	if (loggingVerbosity_ > 0)
		std::cout << "Set Branch Address of susy::Event" << std::endl;
	event = new susy::Event;
	//event->setInput( *inputTree );
	inputTree->SetBranchAddress("susyEvent", &event );

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
	pileupHisto = 0;
}

bool TreeWriter::isData() {
	//event->getEntry(jentry);
	inputTree->GetEntry(1);
	return event->isRealData;
}

void TreeWriter::IncludeAJson(TString const& _fileName) {
	/** Read a Json file which contains good runNumbers and Lumi-sections.
	 * The content will be stored in the class variable 'goodLumiList'.
	 */
	ifstream inputFile(_fileName);
	if(!inputFile.is_open()){
		std::cerr << "Cannot open JSON file " << _fileName << std::endl;
		return;
	}

	std::string line;
	TString jsonText;
	while(inputFile.good()){
		getline(inputFile, line);
		jsonText += line;
	}
	inputFile.close();

	TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
	TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

	TArrayI positions(2);
	positions[1] = 0;
	while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
		TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
		TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

		unsigned run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
		std::set<unsigned>& lumis(goodLumiList[run]);

		TArrayI lumiPos(2);
		lumiPos[1] = 0;
		while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
			TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
			int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
			int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
			for(int lumi(begin); lumi <= end; ++lumi)
				lumis.insert(lumi);
		}
	}
}

void TreeWriter::PileUpWeightFile( std::string const & pileupFileName ) {
	/** Reads the pileup histogram from a given file.
	 */
	TFile *puFile = new TFile( pileupFileName.c_str() );
	pileupHisto = (TH1F*) puFile->Get("pileup");
}

bool TreeWriter::passTrigger() {
	/**
	 * Checks if event passes the HLT trigger paths.
	 *
	 * If the a trigger path contains one of the triggers defined in triggerNames,
	 * true will be returned. If no match is found, false is returned.
	 */
	for( std::vector<const char*>::iterator it = triggerNames.begin();
			it != triggerNames.end(); ++it ) {
		for( susy::TriggerMap::iterator tm = event->hltMap.begin();
				tm != event->hltMap.end(); ++tm ) {
			if ( tm->first.Contains( *it ) && (int(tm->second.second)))
				return true;
		}
	}
	return false;
}

bool TreeWriter::isGoodLumi() const {
	/**
	 * Check if current event is in json file added by 'IncludeAJson(TString )'
	 */
	bool goodLumi = false;
	unsigned run = event->runNumber;
	unsigned lumi = event->luminosityBlockNumber;
	std::map<unsigned, std::set<unsigned> >::const_iterator rItr(goodLumiList.find(run));
	if(rItr != goodLumiList.end()){
		std::set<unsigned>::const_iterator lItr(rItr->second.find(lumi));
		if(lItr != rItr->second.end()) goodLumi = true;
	}
	return goodLumi;
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
		//event->getEntry(jentry);
		inputTree->GetEntry(jentry);

		if ( ! passTrigger() ) continue;
		if ( ! isGoodLumi() ) continue;

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
			if( !(it->isEE() || it->isEB()) && it->momentum.Pt()<20 && it->isEBEtaGap() && it->isEBPhiGap() && it->isEERingGap() && it->isEEDeeGap() && it->isEBEEGap() )
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
			if(!(loose_photon_endcap || loose_photon_barrel || thisphoton.ptJet > 75 ) )
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
		if( photon.size() == 0 ) continue;
		std::sort( photon.begin(), photon.end(), tree::EtGreater);
		if( loggingVerbosity > 1 )
			std::cout << "Found " << photon.size() << " photons" << std::endl;

		// electrons
		std::vector<susy::Electron> eVector = event->electrons["gsfElectrons"];
		for(std::vector<susy::Electron>::iterator it = eVector.begin(); it < eVector.end(); ++it) {
			if( it->momentum.Pt() < 15 || it->momentum.Eta() > 2.6 || !isVetoElectron( *it, *event, loggingVerbosity ) )
				continue;
			tree::Particle thiselectron;
			thiselectron.pt = it->momentum.Pt();
			thiselectron.eta = it->momentum.Eta();
			thiselectron.phi = it->momentum.Phi();
			electron.push_back( thiselectron );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << electron.size() << " electrons" << std::endl;

		// muons
		std::vector<susy::Muon> mVector = event->muons["muons"];
		for( std::vector<susy::Muon>::iterator it = mVector.begin(); it != mVector.end(); ++it) {
			// see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
			if( it->momentum.Pt() < 15 || it->momentum.Eta() > 2.6 || !(it->isPFMuon() && (it->isGlobalMuon() || it->isTrackerMuon())) )
				continue;
			tree::Particle thismuon;
			thismuon.pt = it->momentum.Et();
			thismuon.eta = it->momentum.Eta();
			thismuon.phi = it->momentum.Phi();
			muon.push_back( thismuon );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << muon.size() << " muons" << std::endl;

		// jets
		for(std::vector<susy::PFJet>::iterator it = jetVector.begin();
				it != jetVector.end(); ++it) {

			// compute H_T with uncorrected ak5PFJets, to have larger agreement
			//  with trigger
			if( std::abs( it->momentum.Eta() ) < 3 && it->momentum.Pt() > 40 )
				ht += it->momentum.Pt();

			// scale with JEC
			float scale = 1.;
			if(it->jecScaleFactors.count("L2L3") == 0)
				std::cout << "ERROR: JEC is not available for this jet" << std::endl;
			else
				scale = it->jecScaleFactors.find("L2L3")->second;
			TLorentzVector corrP4 = scale * it->momentum;

			if(std::abs(corrP4.Eta()) > 2.6 ) continue;
			if(corrP4.Pt() < 30 ) continue;
			if( isAdjacentToLepton( *it, electron ) ||  isAdjacentToLepton( *it, muon ) ) continue;
			tree::Jet thisjet;
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
			jet.push_back( thisjet );

			if( loggingVerbosity > 2 )
				std::cout << " p_T, jet = " << thisjet.pt << std::endl;
		}// for jet

		if( jet.size() < 2 )
			continue;
		std::sort( jet.begin(), jet.end(), tree::EtGreater);
		if( loggingVerbosity > 1 )
			std::cout << "Found " << jet.size() << " jets" << std::endl;

		// H_T cut
		if( ht < 450)
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

		// vertices
		nVertex = event->vertices.size();

		tree::Particle thisGenParticle;
		for( std::vector<susy::Particle>::iterator it = event->genParticles.begin(); it != event->genParticles.end(); ++it ) {
			// status 3: particles in matrix element
			// status 2: intermediate particles
			// status 1: final particles (but can decay in geant, etc)
			if( it->momentum.Pt() < 20 || it->status != 1) continue;

			thisGenParticle.pt = it->momentum.Pt();
			thisGenParticle.eta = it->momentum.Eta();
			thisGenParticle.phi = it->momentum.Phi();
			switch( std::abs(it->pdgId) ) {
				case 22: // photon
					genPhoton.push_back( thisGenParticle );
					break;
				case 11: // electron
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

