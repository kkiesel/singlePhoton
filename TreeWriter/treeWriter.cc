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

template <typename Particle>
bool isAdjacentToParticles( const tree::Jet& jet, const std::vector<Particle>& particles, float deltaR_ = 0.3 ) {
	/** Particles near the jet are searched.
	 *
	 * Returns true if a particle is found in a certain radius near the jet.
	 */
	for( typename std::vector<Particle>::const_iterator particle = particles.begin(); particle != particles.end(); ++particle) {
		if ( jet.DeltaR( *particle ) < deltaR_)
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

bool passLooseJetId( const susy::PFJet& jet ) {
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

bool matchLorentzToGenVector( TLorentzVector& lvec, std::vector<tree::Particle>& genParticles, TH2F& hist, float deltaPtRel_ = .3, float deltaR_ = .3 ) {
	bool match = false;
	float dR, dPt;
	TLorentzVector a;
	for( std::vector<tree::Particle>::iterator it = genParticles.begin();
			it != genParticles.end(); ++it ) {
		a.SetPtEtaPhiE( it->pt, it->eta, it->phi, 1  );
		dR = lvec.DeltaR( a );
		dPt = 2*(it->pt-lvec.Pt())/(it->pt+lvec.Pt());
		hist.Fill( dR, dPt );
		if ( dR <= deltaR_ &&  std::abs(dPt) <= deltaPtRel_ )
			match = true;
	}
	return match;
}

void printChildren( int index, susy::ParticleCollection&  particles, int level=0 ){
	/* Prints the daughter and recursivly her's daughters of particle 'index'.
	 * The status of the particle is denoted in paranthesis.
	 */
	for (int i = 0; i< level; ++i )
		std::cout <<"\t";
	std::cout << particles[index].pdgId << " (" << (int)particles[index].status << ")"<< std::endl;

	//susy::ParticleCollection children;
	for( susy::ParticleCollection::iterator it = particles.begin(); it != particles.end(); ++it ) {
		if( it->motherIndex == index )
			printChildren( std::distance(particles.begin(), it ), particles, level+1 );
	}
}

///////////////////////////////////////////////////////////////////////////////
// Here the class implementation begins ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

TreeWriter::TreeWriter( int nFiles, char** fileList, std::string const& outputName) :
	reportEvery(20000),
	processNEvents(-1),
	loggingVerbosity(0),
	splitting(true),
	onlyMetPlots(false),
	hadronicSelection(true),
	photonPtThreshold(80),
	inputTree("susyTree"),
	event(),
	outFile( outputName.c_str(), "recreate" ),
	photonTree("photonTree","Tree for single photon analysis"),
	photonElectronTree("photonElectronTree","Tree for single photon analysis"),
	photonJetTree("photonJetTree","Tree for single photon analysis"),
	eventNumbers("eventNumbers", "Histogram containing number of generated events", 1, 0, 1),
	nPhotons("nPhotons", ";#gamma;#gamma_{jet};#gamma_{e}", 3, -.5, -2.5, 3, -.5, -2.5, 3, -.5, -2.5 )
{

	for( int i = 0; i<nFiles; ++i )
		inputTree.Add( fileList[i] );
	event.setInput( inputTree );

	// Here the number of proceeded events will be stored. For plotting, simply use L*sigma/eventNumber
	eventNumbers.GetXaxis()->SetBinLabel(1,"Number of generated events");

	// Define one dimensional histograms
	hist1D["gMet"] = TH1F("", ";met;", 50, 0, 500 );
	hist1D["eMet"] = TH1F("", ";met;", 50, 0, 500 );
	hist1D["fMet"] = TH1F("", ";met;", 50, 0, 500 );
	hist1D["fMetUp"] = TH1F("", ";met;", 50, 0, 500 );
	hist1D["fMetDown"] = TH1F("", ";met;", 50, 0, 500 );

	// Define two dimensional histograms
	hist2D["matchJet"] = TH2F("", "photon-jet matching;#DeltaR;p_{T, jet}/p_{T, #gamma}", 100, 0, 1, 100, 0, 4 );
	hist2D["matchJetFO"] = TH2F("", "photon-jet matching;#DeltaR;p_{T, jet}/p_{T, #gamma}", 100, 0, 1, 100, 0, 4 );
	hist2D["matchPhoton"] = TH2F("", ";#DeltaR;#Delta p_{T}/p_{T}", 100, 0, 1, 200, -2, 2 );
	hist2D["matchElectron"] = TH2F("", ";#DeltaR;#Delta p_{T}/p_{T}", 100, 0, 1, 200, -2, 2 );
	hist2D["default"] = TH2F("", ";", 1, 0, 1, 1, 0, 1 );

	// Set the keyName as histogram name for one and two dimensional histograms
	for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
			it!= hist1D.end(); ++it ) {
		it->second.SetName( (it->second.GetName() + it->first).c_str() );
		it->second.Sumw2();
	}
	for( std::map<std::string, TH2F>::iterator it = hist2D.begin();
			it!= hist2D.end(); ++it ) {
		it->second.SetName( (it->second.GetName() + it->first).c_str() );
		it->second.Sumw2();
	}
}

TreeWriter::~TreeWriter() {
	/** Deconstructor
	 * Event has to be deleted before the deletion of the tree.
	 */
	event.releaseTree(inputTree);
}

void TreeWriter::SetJsonFile(TString const& filename) {
	/** Read a Json file which contains good runNumbers and Lumi-sections.
	 * The content will be stored in the class variable 'goodLumiList'.
	 */
	if( goodLumiList.size() ) {
		std::cout << "WARNING: Good lumi sections already defined. Now overwritten by " << filename << std::endl;
		return;
	}

	ifstream inputFile( filename );
	if( !inputFile.is_open() )
		std::cerr << "ERROR: Cannot open JSON file " << filename << std::endl;

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

void TreeWriter::SetPileUpWeightFile( std::string const & filename ) {
	/** Reads the pileup histogram from a given file.
	 */

	TFile puFile( filename.c_str() );
	if( puFile.IsZombie() )
		std::cerr << "ERROR: Could not read pileup weight file " << filename << std::endl;

	std::string histogramName = "pileup";
	if( puFile.GetListOfKeys()->Contains( histogramName.c_str() ) )
		pileupHisto = *((TH1F*) puFile.Get( histogramName.c_str() ));
	else
		std::cerr << "ERROR: Could not extract " << histogramName
			<< " histogram from " << filename << std::endl;

	if( loggingVerbosity > 1 )
		std::cout << "Pile-up reweighting histogram added." << std::endl;
}

void TreeWriter::SetQcdWeightFile( std::string const & filename ) {
	/** Reads the pileup histogram from a given file.
	 */
	gSystem->Load("libHistPainter"); // to avoid waring and errors when reading th2 from file

	TFile qcdFile( filename.c_str() );
	if( qcdFile.IsZombie() )
		std::cerr << "ERROR: Could not read pileup weight file " << filename << std::endl;

	std::string histogramName = "qcdWeight";
	if( qcdFile.GetListOfKeys()->Contains( histogramName.c_str() ) )
		qcdWeightHisto = *((TH2F*) qcdFile.Get( histogramName.c_str() ));
	else
		std::cerr << "ERROR: Could not extract " << histogramName
			<< " histogram from " << filename << std::endl;

	if( loggingVerbosity > 1 )
		std::cout << "QCD weighting histogram added." << std::endl;
}

bool TreeWriter::passTrigger() {
	/**
	 * Checks if event passes the HLT trigger paths.
	 *
	 * If the a trigger path contains one of the triggers defined in triggerNames,
	 * true will be returned. If no match is found, false is returned.
	 */
	if( !triggerNames.size() ) {
		std::cout << "WARNING: No triggers selected" << std::endl;
		return true;
	}

	for( std::vector<const char*>::iterator it = triggerNames.begin();
			it != triggerNames.end(); ++it ) {
		for( susy::TriggerMap::iterator tm = event.hltMap.begin();
				tm != event.hltMap.end(); ++tm ) {
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
	if( !goodLumiList.size() ) {
		std::cout << "WARNING: No json file for filtering found" << std::endl;
		return true;
	}
	bool goodLumi = false;
	unsigned run = event.runNumber;
	unsigned lumi = event.luminosityBlockNumber;
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

float TreeWriter::getPileUpWeight(){
	/**
	 * If a pileup weight histogram has been added, the pile-up weight for the
	 * current event is computed.
	 */

	float trueNumInteractions = -1;
	for( susy::PUSummaryInfoCollection::const_iterator iBX = event.pu.begin();
			iBX != event.pu.end() && trueNumInteractions < 0; ++iBX) {
		if (iBX->BX == 0)
			trueNumInteractions = iBX->trueNumInteractions;
	}
	float thisWeight = pileupHisto.GetBinContent( pileupHisto.FindBin( trueNumInteractions ) );

	if( loggingVerbosity > 2 )
		std::cout << "Pile-up weight = " << thisWeight << std::endl;
	return thisWeight;
}

void TreeWriter::getQcdWeights( float pt, float ht_, float & qcdWeight, float & qcdWeightUp, float & qcdWeightDown ){
	if( qcdWeightHisto ) {
		int bin = qcdWeightHisto->FindBin( pt, ht_ );
		float error = qcdWeightHisto->GetBinError(bin);
		qcdWeight = qcdWeightHisto->GetBinContent(bin);
		qcdWeightUp = qcdWeight + error;
		qcdWeightDown = qcdWeight - error;
	} else {
		if( loggingVerbosity > 0 ) std::cout << "WARNING: No qcd weight found." << std::endl;
		qcdWeight = 0;
		qcdWeightUp = 0;
		qcdWeightDown = 0;
	}
}

void TreeWriter::getPtFromMatchedJet( tree::Photon& myPhoton, bool isPhoton=true ) {
	/**
	 * \brief Takes jet p_T as photon p_T
	 *
	 * At first all jets with DeltaR < 0.3 (isolation cone) are searched.
	 * If several jets are found, take the one with the minimal pt difference
	 * compared to the photon. If no such jets are found, keep the photon_pt
	 */
	float returnPt = 0;
	std::vector<short> indices;

	for(std::vector<tree::Jet>::iterator jet = jets.begin();
			jet != jets.end(); ++jet) {

		float deltaR_ = myPhoton.DeltaR( *jet );
		float eRel = jet->pt / myPhoton.pt;
		if( isPhoton )
			hist2D["matchJet"].Fill( deltaR_, eRel );
		else
			hist2D["matchJetFO"].Fill( deltaR_, eRel );

		if (deltaR_ > 0.3 || eRel <= 0.95 ) continue;
		jet->setMatch( tree::kJetAllPhoton );

		if( loggingVerbosity > 2 ) std::cout << " pT_jet / pT_gamma = " << eRel << std::endl;

		indices.push_back( std::distance( jets.begin(), jet ) );
		returnPt = jet->pt;
	}// for jet
	if( indices.size() > 1 ) {
		float pt = 0;
		short minIndex = 0;
		float minPtDifferenz = 1E20; // should be very high

		for( std::vector<short>::const_iterator index = indices.begin();
				index != indices.end(); ++index ) {
			float ptDiff = std::abs(myPhoton.pt - jets.at(*index).pt);
			if ( ptDiff < minPtDifferenz ) {
				minPtDifferenz = ptDiff;
				pt = jets.at(*index).pt;
				minIndex = *index;
			}
		}
		jets.at(minIndex).setMatch( tree::kJetPhoton );
		myPhoton.matchedJetIndex = minIndex;
		myPhoton._ptJet = pt;
	} else if( indices.size() == 1 ) {
		myPhoton.matchedJetIndex = indices.at(0);
		jets.at(0).setMatch( tree::kJetPhoton );
		myPhoton._ptJet = returnPt;
	} else if( loggingVerbosity > 1 )
		std::cout << "No matching jet found, do not change photon_pt." << std::endl;
}

void TreeWriter::fillJets() {
	jets.clear();
	tree::Jet jetToTree;

	std::vector<susy::PFJet> jetVector = event.pfJets.find("ak5")->second;
	for(std::vector<susy::PFJet>::iterator it = jetVector.begin();
			it != jetVector.end(); ++it) {

		TLorentzVector corrP4 = it->jecScaleFactors.at("L1FastL2L3") * it->momentum;
		if( std::abs(corrP4.Eta()) > 3 ) continue;
		if( corrP4.Pt() < 30 ) continue;

		jetToTree.matchInformation = 0;
		if( passLooseJetId( *it ) ) jetToTree.setMatch( tree::kJetGoodId );
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
		jets.push_back( jetToTree );

		if( loggingVerbosity > 2 )
			std::cout << " p_T, jet = " << jetToTree.pt << std::endl;
	}// for jet
	std::sort( jets.begin(), jets.end(), tree::EtGreater);
}

float TreeWriter::getHt() const {
	if( !splitting )
		return 0;

	float returnedHt = 0;
	for(std::vector<tree::Jet>::const_iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {

		if( !jet->isMatch( tree::kJetGoodId ) && !jet->isMatch( tree::kJetPhoton ) ) continue;
		if( jet->pt < 40 || jet->eta > 3. ) continue;

		returnedHt += jet->pt;
	}
	for( std::vector<tree::Photon>::const_iterator photon = photons.begin();
			photon != photons.end(); ++photon ) {
		if( photon->_ptJet == 0 )
			returnedHt += photon->pt;
	}
	for( std::vector<tree::Photon>::const_iterator photon = photonJets.begin();
			photon != photonJets.end(); ++photon ) {
		if( photon->_ptJet == 0 )
			returnedHt += photon->pt;
	}
	for( std::vector<tree::Photon>::const_iterator photon = photonElectrons.begin();
			photon != photonElectrons.end(); ++photon ) {
		if( photon->_ptJet == 0 )
			returnedHt += photon->pt;
	}

	return returnedHt;
}

unsigned int TreeWriter::countGoodJets( bool clean ) const {
	unsigned int number = 0;
	for(std::vector<tree::Jet>::const_iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {
		if( jet->pt < 30 || jet->eta > 2.5 ) continue;
		if( !jet->isMatch( tree::kJetGoodId ) ) continue;

		if( isAdjacentToParticles<tree::Particle>( *jet, electrons ) ) continue;
		if( isAdjacentToParticles<tree::Particle>( *jet, muons ) ) continue;
		if( clean ) {
			if( isAdjacentToParticles<tree::Photon>( *jet, photons ) ) continue;
			if( isAdjacentToParticles<tree::Photon>( *jet, photonElectrons ) ) continue;
			if( isAdjacentToParticles<tree::Photon>( *jet, photonJets ) ) continue;
		}
		++number;
	}
	return number;
}
void TreeWriter::SetBranches( TTree& tree ) {
	tree.Branch("jets", &jets);
	tree.Branch("electrons", &electrons);
	tree.Branch("muons", &muons);
	tree.Branch("genPhotons", &genPhotons);
	tree.Branch("genElectrons", &genElectrons);
	tree.Branch("met", &met, "met/F");
	tree.Branch("type0met", &type0met, "type0met/F");
	tree.Branch("type1met", &type1met, "type1met/F");
	tree.Branch("ht", &ht, "ht/F");
	tree.Branch("weight", &weight, "weight/F");
	tree.Branch("nVertex", &nVertex, "nVertex/I");
	tree.Branch("nGoodJets", &nGoodJets, "nGoodJets/i");
	tree.Branch("runNumber", &runNumber, "runNumber/i");
	tree.Branch("eventNumber", &eventNumber, "eventNumber/i");
	tree.Branch("luminosityBlockNumber", &luminosityBlockNumber, "luminosityBlockNumber/i");
}

void TreeWriter::Loop() {
	/**
	 * \brief Loops over input chain and fills tree
	 *
	 * This is the major function of treeWriter, which initialize the output, loops
	 * over all events and fill the tree. In the end, the tree is saved to the
	 * output File
	 */

	// get number of events to be proceeded
	Long64_t nentries = inputTree.GetEntries();
	// store them in histo
	eventNumbers.Fill( "Number of generated events", nentries );

	if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;
	if( loggingVerbosity > 0 )
		std::cout << "Processing " << processNEvents << " ouf of "
			<< nentries << " events. " << std::endl;

	photonTree.Branch( "photons", &photons );
	photonElectronTree.Branch( "photons", &photonElectrons );
	photonJetTree.Branch( "photons", &photonJets );
	SetBranches( photonTree );
	SetBranches( photonElectronTree );
	SetBranches( photonJetTree );

	// Declaration for objects saved in Tree
	tree::Photon photonToTree;
	tree::Particle electronToTree;
	tree::Particle muonToTree;

	for (long jentry=0; jentry < processNEvents; ++jentry) {
		if ( loggingVerbosity>1 || jentry%reportEvery==0 ) std::cout << jentry << " / " << processNEvents << std::endl;
		event.getEntry(jentry);

		bool printCascade = false;
		for( susy::ParticleCollection::iterator it = event->genParticles.begin(); printCascade && it != event->genParticles.end(); ++it ){
			if( it->motherIndex == -1 ){
				printChildren( std::distance(event.genParticles.begin(), it ), event.genParticles );
			}
		}

		if ( event.isRealData )
			if ( !isGoodLumi() || !passTrigger() || !event.passMetFilters() ) continue;

		// vertices
		nVertex = numberOfGoodVertexInCollection( event.vertices );
		if( !nVertex ) continue;

		photons.clear();
		photonJets.clear();
		photonElectrons.clear();
		jets.clear();
		electrons.clear();
		muons.clear();
		genElectrons.clear();
		genPhotons.clear();
		runNumber = event.runNumber;
		eventNumber = event.eventNumber;
		luminosityBlockNumber = event.luminosityBlockNumber;
		weight = event.isRealData ? 1. : getPileUpWeight();

		// genParticles
		tree::Particle thisGenParticle;
		for( std::vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); ++it ) {

			// status 3: particles in matrix element
			// status 2: intermediate particles
			// status 1: final particles (but can decay in geant, etc)
			// TODO: Check substraction of 10
			if( it->momentum.Pt() < photonPtThreshold-10 || it->status != 1) continue;

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

		// electrons
		std::vector<susy::Electron> eVector = event.electrons["gsfElectrons"];
		for(std::vector<susy::Electron>::iterator it = eVector.begin(); it < eVector.end(); ++it) {
			if( it->momentum.Pt() < 15 || std::abs(it->momentum.Eta()) > 2.6 || !isVetoElectron( *it, event, loggingVerbosity ) )
				continue;
			electronToTree.pt = it->momentum.Pt();
			electronToTree.eta = it->momentum.Eta();
			electronToTree.phi = it->momentum.Phi();
			electrons.push_back( electronToTree );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << electrons.size() << " electrons" << std::endl;

		// muons
		std::vector<susy::Muon> mVector = event.muons["muons"];
		for( std::vector<susy::Muon>::iterator it = mVector.begin(); it != mVector.end(); ++it) {
			// see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
			if( it->momentum.Pt() < 15 || std::abs(it->momentum.Eta()) > 2.6 || !(it->isPFMuon() && (it->isGlobalMuon() || it->isTrackerMuon())) )
				continue;
			muonToTree.pt = it->momentum.Et();
			muonToTree.eta = it->momentum.Eta();
			muonToTree.phi = it->momentum.Phi();
			muons.push_back( muonToTree );
		}
		if( loggingVerbosity > 1 )
			std::cout << "Found " << muons.size() << " muons" << std::endl;

		// get jets
		fillJets();

		// photons
		std::vector<susy::Photon> photonVector = event.photons["photons"];
		for(std::vector<susy::Photon>::iterator it = photonVector.begin();
				it != photonVector.end(); ++it ) {
			float eta = std::abs( it->momentum.Eta() );
			if( it->momentum.Pt() < photonPtThreshold || eta >= susy::etaGapBegin )
				continue;
			photonToTree.chargedIso = chargedHadronIso_corrected(*it, event.rho25);
			photonToTree.neutralIso = neutralHadronIso_corrected(*it, event.rho25);
			photonToTree.photonIso = photonIso_corrected(*it, event.rho25);
			photonToTree.pt = it->momentum.Pt();
			photonToTree.eta = it->momentum.Eta();
			photonToTree.phi = it->momentum.Phi();
			photonToTree.r9 = it->r9;
			photonToTree.sigmaIetaIeta = it->sigmaIetaIeta;
			photonToTree.sigmaIphiIphi = it->sigmaIphiIphi;
			photonToTree.hadTowOverEm = it->hadTowOverEm;
			photonToTree.pixelseed = it->nPixelSeeds;
			photonToTree.conversionSafeVeto = it->passelectronveto;
			photonToTree.genInformation = 0;
			if( matchLorentzToGenVector( it->momentum, genPhotons, hist2D["matchPhoton"], 1e6, .05 ) )
				photonToTree.setGen( tree::kGenPhoton );
			if( matchLorentzToGenVector( it->momentum, genElectrons, hist2D["matchElectron"], 1e6, .05 ) )
				photonToTree.setGen( tree::kGenElectron );

			if( matchLorentzToGenVector( it->momentum, electrons, hist2D["default"], 1e6 ) ||
					matchLorentzToGenVector( it->momentum, muons, hist2D["default"], 1e6 ) )
				photonToTree.setGen( tree::kNearLepton );


			getPtFromMatchedJet( photonToTree, photonToTree.isGen( tree::kGenPhoton) );

			if( splitting ) {

				//photon definition barrel
				bool isPhotonOrElectron = eta < susy::etaGapBegin
					&& photonToTree.hadTowOverEm < 0.05
					&& photonToTree.sigmaIetaIeta < 0.012
					&& photonToTree.chargedIso < 2.6
					&& photonToTree.neutralIso < 3.5+0.04*photonToTree.pt
					&& photonToTree.photonIso < 1.3+0.005*photonToTree.pt;

				bool isPhotonJet = eta < susy::etaGapBegin
					&& !photonToTree.pixelseed
					&& photonToTree.hadTowOverEm < 0.05
					&& photonToTree.sigmaIetaIeta < 0.012
					&& photonToTree.chargedIso < 26 && photonToTree.chargedIso > 0.26
					&& photonToTree.neutralIso < 35+0.4*photonToTree.pt && photonToTree.neutralIso > 0.35+0.004*photonToTree.neutralIso
					&& photonToTree.photonIso < 13+0.05*photonToTree.pt && photonToTree.photonIso > 0.13+0.0005*photonToTree.pt;

					if( isPhotonOrElectron ) {
						if( photonToTree.pixelseed )
							photonElectrons.push_back( photonToTree );
						else
							photons.push_back( photonToTree );
					} else if ( isPhotonJet )
						photonJets.push_back( photonToTree );

			} else // no splitting, put everything in the vector 'photons'
				photons.push_back( photonToTree );

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
		nPhotons.Fill( photons.size(), photonJets.size(), photonElectrons.size() );

		// met
		met = event.metMap["pfMet"].met();
		type0met = event.metMap["pfType01CorrectedMet"].met();
		type1met = event.metMap["pfType1CorrectedMet"].met();
		if( loggingVerbosity > 2 )
			std::cout << " met = " << met << std::endl;

		ht = getHt();
		nGoodJets = countGoodJets( splitting );

		if( splitting && hadronicSelection && ( nGoodJets < 2 || ht < 500 ) ) continue;

		if( splitting ) {
			float gPt = photons.size()         ? photons.at(0).pt         : 0;
			float ePt = photonElectrons.size() ? photonElectrons.at(0).pt : 0;
			float fPt = photonJets.size()      ? photonJets.at(0).pt      : 0;

			bool isPhotonEvent = false;
			bool isPhotonJetEvent = false;
			bool isPhotonElectronEvent = false;

			if( gPt && gPt > fPt && gPt > ePt )
				isPhotonEvent = true;
			if( ePt && ePt > fPt && ePt > gPt )
				isPhotonElectronEvent = true;
			if( fPt && fPt > ePt && fPt > gPt )
				isPhotonJetEvent = true;

			if( isPhotonEvent ) {
				photonTree.Fill();
				hist1D["gMet"].Fill( met, weight );
			}
			if( isPhotonJetEvent) {
				photonJetTree.Fill();
				float qcdWeight=0, qcdWeightUp=0, qcdWeightDown=0;
				getQcdWeights( photonJets.at(0).ptJet(), ht, qcdWeight, qcdWeightUp, qcdWeightDown );
				hist1D["fMet"].Fill( met, weight*qcdWeight );
				hist1D["fMetUp"].Fill( met, weight*qcdWeightUp );
				hist1D["fMetDown"].Fill( met, weight*qcdWeightDown );
			}
			if( isPhotonElectronEvent ) {
				photonElectronTree.Fill();
				hist1D["eMet"].Fill( met, weight );
			}
			if( isPhotonEvent+isPhotonElectronEvent+isPhotonJetEvent > 1 )
				std::cout <<"ERROR: One event is control and signal at once!" << std::endl;
		} else if( photons.size() ) // no splitting
				photonTree.Fill();

	} // for jentry

	outFile.cd();
	if( !onlyMetPlots ) {
		photonTree.Write();
		if( splitting ) {
			photonElectronTree.Write();
			photonJetTree.Write();
			nPhotons.Write();
		}
		eventNumbers.Write();
		for( std::map<std::string, TH2F>::iterator it = hist2D.begin();
				it!= hist2D.end(); ++it )
			it->second.Write();
	}
	TPRegexp expFilename( ".*/tree_([0-9]+_[0-9]+)_375.root" );
	TObjArray *arr = expFilename.MatchS( inputTree.GetCurrentFile()->GetName() );
	std::string histoNameAppendix = "";
	if( arr->GetLast() >0 )
		histoNameAppendix = (std::string)(((TObjString *)arr->At(1))->GetString());
	else
		std::cout << "Could not extract grid parameters from filename." << std::endl;

	for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
			it!= hist1D.end(); ++it ) {
		it->second.SetName( (it->second.GetName() + histoNameAppendix ).c_str() );
		it->second.Write();
	}
}

