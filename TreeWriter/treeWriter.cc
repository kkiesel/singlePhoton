#include "treeWriter.h"

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

template <typename Particle>
bool isAdjacentToParticles( const susy::PFJet& jet, const std::vector<Particle>& particles, float deltaR_ = 0.3 ) {
	/** Particles near the jet are searched.
	 *
	 * Returns true if a particle is found in a certain radius near the jet.
	 */
	TLorentzVector a;
	for( typename std::vector<Particle>::const_iterator particle = particles.begin(); particle != particles.end(); ++particle) {
		a.SetPtEtaPhiE( 1,particle->eta, particle->phi,1  );
		if ( jet.momentum.DeltaR( a ) < deltaR_)
			return true;
	}
	return false;
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
	susy::Track track = event.tracks[electron.gsfTrackIndex];
	float d0 = track.d0();
	float dZ = track.vertex.Z();
	bool isElectron = false;
	float eta = std::abs(electron.momentum.Eta());
	isElectron  = (
		eta < susy::etaGapBegin
			&& ( fabs(electron.deltaEtaSuperClusterTrackAtVtx) < 0.007
				|| fabs(electron.deltaPhiSuperClusterTrackAtVtx) < 0.8
				|| electron.sigmaIetaIeta < 0.01
				|| electron.hcalOverEcalBc < 0.15
				|| d0 < 0.04
				|| dZ < 0.2
				|| iso < 0.15 )
		)||( ( susy::etaGapEnd < eta && eta < susy::etaMax )
			&& ( fabs(electron.deltaEtaSuperClusterTrackAtVtx) < 0.01
				|| fabs(electron.deltaPhiSuperClusterTrackAtVtx) < 0.7
				|| electron.sigmaIetaIeta < 0.03
				|| d0 < 0.04
				|| dZ < 0.2
				|| iso < 0.15 )
		);
	return isElectron;
}

bool looseJetId( const susy::PFJet& jet ) {
	/**
	 * \brief Apply loose cut on jets.
	 *
	 * See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_7_TeV_data_a
	 * for more information.
	 */
	double energy = jet.momentum.E();
	return jet.neutralHadronEnergy / energy < 0.99
			&& jet.neutralEmEnergy / energy < 0.99
			&& jet.nConstituents > 1
			&& ( std::abs(jet.momentum.Eta())>=2.4
				|| ( jet.chargedHadronEnergy / energy > 0
					&& jet.chargedMultiplicity > 0
					&& jet.chargedEmEnergy < 0.99 ) );
}

bool goodVertex( susy::Vertex& vtx ) {
	/** Definition of a good vertex
	 */
	return (!vtx.isFake() &&
		vtx.ndof > 4 &&
		std::abs((vtx.position).z()) < 24.0 &&
		std::abs((vtx.position).Perp()) < 2.0 );
}

unsigned int numberOfGoodVertexInCollection( std::vector<susy::Vertex>& vertexVector ) {
	unsigned int number = 0;
	for( std::vector<susy::Vertex>::iterator vtx = vertexVector.begin();
			vtx != vertexVector.end(); ++vtx ) {
		if( goodVertex( *vtx ) )
			number++;
	}
	return number;
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
	delete photonTree;
	delete photonElectronTree;
	delete photonJetTree;
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
	event->setInput( *inputTree );

	// Here the number of proceeded events will be stored. For plotting, simply use L*sigma/eventNumber
	eventNumbers = new TH1F("eventNumbers", "Histogram containing number of generated events", 1, 0, 1);
	eventNumbers->GetXaxis()->SetBinLabel(1,"Number of generated events");
	matchingHisto = new TH2F("matchingPhotonJet", "photon-jet matching;#DeltaR;p_{T, jet}/p_{T, #gamma}", 1000, 0, 3, 1000, 0, 4 );

	// open the output file
	if (loggingVerbosity_>0)
		std::cout << "Open file " << outputName << " for writing." << std::endl;
	outFile = new TFile( outputName.c_str(), "recreate" );
	photonTree = new TTree("photonTree","Tree for single photon analysis");
	photonElectronTree = new TTree("photonElectronTree","Tree for single photon analysis");
	photonJetTree = new TTree("photonJetTree","Tree for single photon analysis");

	// set default parameter
	processNEvents = -1;
	reportEvery = 1000;
	loggingVerbosity = loggingVerbosity_;
	pileupHisto = 0;
	useAdditionalFOCut = false;
}

bool TreeWriter::isData() {
	//event->getEntry(jentry);
	event->getEntry(1);
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
	if( loggingVerbosity > 1 )
		std::cout << "JSON file for filtering included." << std::endl;
}

void TreeWriter::PileUpWeightFile( std::string const & pileupFileName ) {
	/** Reads the pileup histogram from a given file.
	 */
	TFile *puFile = new TFile( pileupFileName.c_str() );
	pileupHisto = (TH1F*) puFile->Get("pileup");
	if( loggingVerbosity > 1 )
		std::cout << "Pile-up reweighting histogram added." << std::endl;
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
			if ( tm->first.Contains( *it ) && (int(tm->second.second))) {
				return true;
				if( loggingVerbosity > 1 )
					std::cout << "Pass trigger requirement." << std::endl;
			}
		}
	}
	if( loggingVerbosity > 1 )
		std::cout << "Fail trigger requirement." << std::endl;
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
	if( loggingVerbosity > 1 && goodLumi )
		std::cout << "Event is in a good lumi section." << std::endl;
	if( loggingVerbosity > 1 && !goodLumi )
		std::cout << "Event is not in a good lumi section." << std::endl;
	return goodLumi;
}

float TreeWriter::getPileUpWeight() const {
	float thisWeight = 1;
	if (pileupHisto == 0) {
		thisWeight = 1.;
	} else {
		float trueNumInteractions = -1;
		for( susy::PUSummaryInfoCollection::const_iterator iBX = event->pu.begin();
				iBX != event->pu.end() && trueNumInteractions < 0; ++iBX) {
			if (iBX->BX == 0)
				trueNumInteractions = iBX->trueNumInteractions;
		}
		thisWeight = pileupHisto->GetBinContent( pileupHisto->FindBin( trueNumInteractions ) );
	}
	if( loggingVerbosity > 2 )
		std::cout << "Pile-up weight = " << thisWeight << std::endl;
	return thisWeight;
}

float TreeWriter::getPtFromMatchedJet( const susy::Photon& myPhoton, bool fillHisto=true ) const {
	/**
	 * \brief Takes jet p_T as photon p_T
	 *
	 * At first all jets with DeltaR < 0.3 (isolation cone) are searched.
	 * If several jets are found, take the one with the minimal pt difference
	 * compared to the photon. If no such jets are found, keep the photon_pt
	 */
	std::vector<susy::PFJet> jetColl = event->pfJets["ak5"];
	std::vector<susy::PFJet> nearJets;

	for(std::vector<susy::PFJet>::const_iterator it = jetColl.begin();
			it != jetColl.end(); ++it) {
		float scale = 1.;
		std::map<TString,Float_t>::const_iterator s_it = it->jecScaleFactors.find("L1FastL2L3");
		if (s_it == it->jecScaleFactors.end()) {
			std::cout << "JEC is not available for this jet!!!" << std::endl;
		} else {
			scale = s_it->second;
		}
		TLorentzVector corrP4 = scale * it->momentum;
		float deltaR_ = myPhoton.momentum.DeltaR( corrP4 );
		float eRel = corrP4.Pt() / myPhoton.momentum.Pt();
		if( fillHisto )
			matchingHisto->Fill( deltaR_, eRel );
		if (deltaR_ > 0.3 || eRel <= 0.95 ) continue;
		if( loggingVerbosity > 2 )
			std::cout << " pT_jet / pT_gamma = " << eRel << std::endl;
		nearJets.push_back( *it );
	}// for jet

	if ( nearJets.size() == 0 ) {
		if( loggingVerbosity > 1 )
			std::cout << "No matching jet found, do not change photon_pt." << std::endl;
		return 0;
	} else if ( nearJets.size() == 1 )
		return nearJets.at(0).momentum.Pt();
	else {
		std::cout << "More than one jet found, set the photon_ptJet to the nearest value in pt." << std::endl;
		float pt = 0;
		float minPtDifferenz = 1E20; // should be very high
		for( std::vector<susy::PFJet>::iterator it = nearJets.begin(), jetEnd = nearJets.end();
				it != jetEnd; ++it ) {
			float ptDiff = std::abs(myPhoton.momentum.Pt() - it->momentum.Pt());
			if (  ptDiff < minPtDifferenz ) {
				minPtDifferenz = ptDiff;
				pt = it->momentum.Pt();
			}
		}
		return pt;
	}
}

std::vector<tree::Jet> TreeWriter::getJets( const std::vector<tree::Photon>& excludedPhotons ) const {
	tree::Jet jetToTree;
	std::vector<tree::Jet> returnedJets;

	std::vector<susy::PFJet> jetVector = event->pfJets["ak5"];
	for(std::vector<susy::PFJet>::iterator it = jetVector.begin();
			it != jetVector.end(); ++it) {
		if( !looseJetId( *it ) ) continue;
		if( !it->passPuJetIdLoose( susy::kPUJetIdFull ) ) continue; // has to be changed for sure

		// scale with JEC
		float scale = 1.;
		if(it->jecScaleFactors.count("L1FastL2L3") == 0)
			std::cout << "ERROR: JEC is not available for this jet" << std::endl;
		else
			scale = it->jecScaleFactors.find("L1FastL2L3")->second;
		TLorentzVector corrP4 = scale * it->momentum;

		if( std::abs(corrP4.Eta()) > 2.6 ) continue;
		if( corrP4.Pt() < 30 ) continue;
		if( isAdjacentToParticles<tree::Particle>( *it, electrons ) ) continue;
		if( isAdjacentToParticles<tree::Particle>( *it, muons ) ) continue;
		if( isAdjacentToParticles<tree::Photon>( *it, excludedPhotons ) ) continue;

		jetToTree.pt = corrP4.Pt();
		jetToTree.eta = corrP4.Eta();
		jetToTree.phi = corrP4.Phi();
		jetToTree.bCSV = it->bTagDiscriminators[susy::kCSV];
		// jet composition
		jetToTree.chargedHadronEnergy = it->chargedHadronEnergy;
		jetToTree.neutralHadronEnergy = it->neutralHadronEnergy;
		jetToTree.photonEnergy = it->photonEnergy;
		jetToTree.electronEnergy = it->electronEnergy;
		jetToTree.muonEnergy = it->muonEnergy;
		jetToTree.HFHadronEnergy = it->HFHadronEnergy;
		jetToTree.HFEMEnergy = it->HFEMEnergy;
		jetToTree.chargedEmEnergy = it->chargedEmEnergy;
		jetToTree.chargedMuEnergy = it->chargedMuEnergy;
		jetToTree.neutralEmEnergy = it->neutralEmEnergy;
		returnedJets.push_back( jetToTree );

		if( loggingVerbosity > 2 )
			std::cout << " p_T, jet = " << jetToTree.pt << std::endl;
	}// for jet
	std::sort( returnedJets.begin(), returnedJets.end(), tree::EtGreater);
	return returnedJets;
}

float TreeWriter::getHt( const tree::Photon& photon ) const {
	// ht
	float returnedHt = 0;
	std::vector<susy::PFJet> jetVector = event->pfJets["ak5"];
	for(std::vector<susy::PFJet>::iterator it = jetVector.begin();
			it != jetVector.end(); ++it) {

		if( !looseJetId( *it ) ) continue;
		if( !it->passPuJetIdLoose( 0 ) ) continue; // has to be changed for sure

		// scale with JEC
		float scale = 1.;
		if(it->jecScaleFactors.count("L1FastL2L3") == 0)
			std::cout << "ERROR: JEC is not available for this jet" << std::endl;
		else
			scale = it->jecScaleFactors.find("L1FastL2L3")->second;
		TLorentzVector corrP4 = scale * it->momentum;

		if( corrP4.Pt() < 40 || corrP4.Eta() > 3. )
			continue;

		returnedHt += corrP4.Pt();
		if( photon._ptJet == 0 )
			returnedHt += photon.pt;
	}
	return returnedHt;
}

void TreeWriter::SetBranches( TTree& tree ) {
	tree.Branch("jets", &jets);
	tree.Branch("electrons", &electrons);
	tree.Branch("muons", &muons);
	tree.Branch("met", &met, "met/F");
	tree.Branch("metPhi", &met_phi, "metPhi/F");
	tree.Branch("type1met", &type1met, "type1met/F");
	tree.Branch("type1metPhi", &type1met_phi, "type1metPhi/F");
	tree.Branch("ht", &ht, "ht/F");
	tree.Branch("nVertex", &nVertex, "nVertex/I");
	tree.Branch("weight", &weight, "weight/D");
	tree.Branch("runNumber", &runNumber, "runNumber/i");
	tree.Branch("eventNumber", &eventNumber, "eventNumber/i");
	tree.Branch("luminosityBlockNumber", &luminosityBlockNumber, "luminosityBlockNumber/i");
	tree.Branch("genElectrons", &genElectrons);
	tree.Branch("genPhotons", &genPhotons);
}

bool TreeWriter::passRecommendedMetFilters() const {
	bool pass =  event->passMetFilter( susy::kCSCBeamHalo )
		&& event->passMetFilter( susy::kHcalNoise )
		&& event->passMetFilter( susy::kHcalLaserOccupancy )
		&& event->passMetFilter( susy::kEcalDeadCellTP )
		&& event->passMetFilter( susy::kTrackingFailure )
		&& event->passMetFilter( susy::kEEBadSC )
		&& event->passMetFilter( susy::kEcalLaserCorr ) // optional for 22Jan2013
		// Tracking odd events filters (tracking POG filters)
		&& event->passMetFilter( susy::kManyStripClus53X )
		&& event->passMetFilter( susy::kTooManyStripClus53X )
		&& event->passMetFilter( susy::kLogErrorTooManyClusters );
	if( loggingVerbosity > 1 && pass )
		std::cout << "Passed MET filters." << std::endl;
	if( loggingVerbosity > 1 && !pass)
		std::cout << "Failed MET filters." << std::endl;
	return pass;
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

	photonTree->Branch( "photons", &photons );
	photonElectronTree->Branch( "photonElectrons", &photonElectrons );
	photonJetTree->Branch( "photonJets", &photonJets );
	SetBranches( *photonTree );
	SetBranches( *photonElectronTree );
	SetBranches( *photonJetTree );

	// Declaration for objects saved in Tree
	tree::Photon photonToTree;
	tree::Particle electronToTree;
	tree::Particle muonToTree;

	for (long jentry=0; jentry < processNEvents; ++jentry) {
		if ( loggingVerbosity>1 || jentry%reportEvery==0 ) std::cout << jentry << " / " << processNEvents << std::endl;
		event->getEntry(jentry);

		if ( event->isRealData )
			if ( !isGoodLumi() || !passTrigger() ) continue;
		if( !passRecommendedMetFilters() ) continue;

		// vertices
		nVertex = numberOfGoodVertexInCollection( event->vertices );
		if( !nVertex ) continue;

		photons.clear();
		photonJets.clear();
		photonElectrons.clear();
		jets.clear();
		electrons.clear();
		muons.clear();
		genElectrons.clear();
		genPhotons.clear();
		ht = 0;
		runNumber = event->runNumber;
		eventNumber = event->eventNumber;
		luminosityBlockNumber = event->luminosityBlockNumber;
		weight = getPileUpWeight();

		// photons
		std::vector<susy::Photon> photonVector = event->photons["photons"];

		for(std::vector<susy::Photon>::iterator it = photonVector.begin();
				it != photonVector.end(); ++it ) {
			float eta = std::abs( it->momentum.Eta() );
			if( it->momentum.Pt() < 80 || eta >= susy::etaGapBegin )
				continue;
			photonToTree.chargedIso = chargedHadronIso_corrected(*it, event->rho25);
			photonToTree.neutralIso = neutralHadronIso_corrected(*it, event->rho25);
			photonToTree.photonIso = photonIso_corrected(*it, event->rho25);
			photonToTree.pt = it->momentum.Pt();
			photonToTree.eta = it->momentum.Eta();
			photonToTree.phi = it->momentum.Phi();
			photonToTree.r9 = it->r9;
			photonToTree.sigmaIetaIeta = it->sigmaIetaIeta;
			photonToTree.hadTowOverEm = it->hadTowOverEm;
			photonToTree.pixelseed = it->nPixelSeeds;
			photonToTree.conversionSafeVeto = it->passelectronveto;
			photonToTree.genInformation = 0;

			//photon definition barrel
			bool isPhotonOrElectron = eta < susy::etaGapBegin
				&& it->hadTowOverEm<0.05
				&& it->sigmaIetaIeta<0.012
				&& photonToTree.chargedIso<2.6
				&& photonToTree.neutralIso<3.5+0.04*photonToTree.pt
				&& photonToTree.photonIso<1.3+0.005*photonToTree.pt;

			bool additionalFOCut = eta < susy::etaGapBegin
				&& it->hadTowOverEm<0.05
				&& it->sigmaIetaIeta<0.014
				&& photonToTree.chargedIso<14
				&& photonToTree.neutralIso<3.5+0.04*photonToTree.pt
				&& photonToTree.photonIso<1.3+0.005*photonToTree.pt
				&& ( it->sigmaIetaIeta >=0.012 || photonToTree.chargedIso>=2.6 );

			if( isPhotonOrElectron ) {
				photonToTree._ptJet = getPtFromMatchedJet( *it );
				if( photonToTree.pixelseed ) {
					photonElectrons.push_back( photonToTree );
				} else {
					photons.push_back( photonToTree );
				}
			} else if ( !useAdditionalFOCut || additionalFOCut ) {
				photonToTree._ptJet = getPtFromMatchedJet( *it );
				photonJets.push_back( photonToTree );
			}

			if( loggingVerbosity > 2 )
				std::cout << " p_T, gamma = " << photonToTree.pt << std::endl;
		}
		std::sort( photons.begin(), photons.end(), tree::EtGreater );
		std::sort( photonElectrons.begin(), photonElectrons.end(), tree::EtGreater );
		std::sort( photonJets.begin(), photonJets.end(), tree::EtGreater );
		if( loggingVerbosity > 1 )
			std::cout << "Found " << photons.size() << " photons, "
					<< photonJets.size() << " photon_{jets} and "
					<< photonElectrons.size() << " photon electrons." << std::endl;

		// electrons
		std::vector<susy::Electron> eVector = event->electrons["gsfElectrons"];
		for(std::vector<susy::Electron>::iterator it = eVector.begin(); it < eVector.end(); ++it) {
			if( it->momentum.Pt() < 15 || it->momentum.Eta() > 2.6 || !isVetoElectron( *it, *event, loggingVerbosity ) )
				continue;
			electronToTree.pt = it->momentum.Pt();
			electronToTree.eta = it->momentum.Eta();
			electronToTree.phi = it->momentum.Phi();
			electrons.push_back( electronToTree );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << electrons.size() << " electrons" << std::endl;

		// muons
		std::vector<susy::Muon> mVector = event->muons["muons"];
		for( std::vector<susy::Muon>::iterator it = mVector.begin(); it != mVector.end(); ++it) {
			// see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
			if( it->momentum.Pt() < 15 || it->momentum.Eta() > 2.6 || !(it->isPFMuon() && (it->isGlobalMuon() || it->isTrackerMuon())) )
				continue;
			muonToTree.pt = it->momentum.Et();
			muonToTree.eta = it->momentum.Eta();
			muonToTree.phi = it->momentum.Phi();
			muons.push_back( muonToTree );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << muons.size() << " muons" << std::endl;

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
					genPhotons.push_back( thisGenParticle );
					break;
				case 11: // electron
					genElectrons.push_back( thisGenParticle );
					break;
			}
		}

		if( photons.size() ) {
			ht = getHt( photons.at(0) );
			jets = getJets( photons );
			if( ht < 450 || jets.size() < 2 ) continue;
			photonTree->Fill();
		} else if( photonElectrons.size() ) {
			ht = getHt( photonElectrons.at(0) );
			jets = getJets( photonElectrons );
			if( ht < 450 || jets.size() < 2 ) continue;
			photonElectronTree->Fill();
		} else if( photonJets.size() ) {
			ht = getHt( photonJets.at(0) );
			jets = getJets( photonJets );
			if( ht < 450 || jets.size() < 2) continue;
			photonJetTree->Fill();
		}

	} // for jentry

	outFile->cd();
	eventNumbers->Write();
	matchingHisto->Write();
	photonTree->Write();
	photonElectronTree->Write();
	photonJetTree->Write();
	outFile->Write();
	outFile->Close();
}

