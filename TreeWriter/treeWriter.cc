#include "treeWriter.h"

using namespace std;

float effectiveAreaElectron( float eta ) {
	/** Returns the effective area for the isolation criteria for electrons.
	 * See https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection
	 * only for Delta R = 0.4 on 2012 Data
	 */
	eta = fabs( eta );
	float ea;

	if( eta < 1.0 ) ea = 0.208;
	else if( eta < 1.479 ) ea = 0.209;
	else if( eta < 2.0 ) ea = 0.115;
	else if( eta < 2.2 ) ea = 0.143;
	else if( eta < 2.3 ) ea = 0.183;
	else if( eta < 2.4 ) ea = 0.194;
	else ea = 0.261;

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

bool goodVertex( const susy::Vertex& vtx ) {
	/** Definition of a good vertex. Returns true if the vertex is good.
	 */
	return (!vtx.isFake() &&
		vtx.ndof > 4 &&
		std::abs((vtx.position).z()) < 24.0 &&
		std::abs((vtx.position).Perp()) < 2.0 );
}

unsigned int numberOfGoodVertexInCollection( const std::vector<susy::Vertex>& vertexVector ) {
	/* Counts the number of good vertices in the vertex Vector
	 */
	unsigned int number = 0;
	for( std::vector<susy::Vertex>::const_iterator vtx = vertexVector.begin();
			vtx != vertexVector.end(); ++vtx ) {
		if( goodVertex( *vtx ) )
			number++;
	}
	return number;
}

bool matchLorentzToGenVector( TLorentzVector& lvec, std::vector<tree::Particle>& genParticles, TH2F* hist=NULL, float deltaR_=.3, float deltaPtRel_=1e6 ) {
	/* Try to match a TLorentzVector any element in a given vector.
	 *
	 * lvec: vector for which match is computed
	 * genParticles: vector of many particles for which a match to lvec is checked.
	 * hist: fill the histogram with dR and relPt
	 */
	bool match = false;
	float dR, dPt;
	TLorentzVector a;
	for( std::vector<tree::Particle>::const_iterator it = genParticles.begin();
			it != genParticles.end(); ++it ) {
		a.SetPtEtaPhiE( it->pt, it->eta, it->phi, 1  );
		dR = lvec.DeltaR( a );
		dPt = 2*(it->pt-lvec.Pt())/(it->pt+lvec.Pt());

		if( hist ) hist->Fill( dR, dPt );

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
	std::cout << particles[index].pdgId << " (" << (int)particles[index].status << ") " << particles[index].momentum.Pt() << std::endl;

	//susy::ParticleCollection children;
	for( susy::ParticleCollection::const_iterator it = particles.begin(); it != particles.end(); ++it ) {
		if( it->motherIndex == index )
			printChildren( std::distance<susy::ParticleCollection::const_iterator>(particles.begin(), it ), particles, level+1 );
	}
}

void fillMetFilterBitHistogram( TH1F& hist, int filterBit ) {
	/* Store the information about the filtered events in a histogram.
	 *
	 * If the event is clean, the 0 bin will be incremented.
	 * The other filters are filled in the order in susy::MetFilters + 1
	 */
	bool eventClean = true;
	for( int filterIndex = 0; filterIndex < susy::nMetFilters; ++filterIndex )
		if( !(filterBit & (1 << filterIndex)) ) {
			hist.AddBinContent( filterIndex+1 );
			eventClean = false;
		}
	if( eventClean )
		hist.AddBinContent( 0 );
}

///////////////////////////////////////////////////////////////////////////////
// Here the class implementation begins ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

TreeWriter::TreeWriter( int nFiles, char** fileList, std::string const& outputName) :
	reportEvery(200000),
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
	nPhotons("nPhotons", ";#gamma;#gamma_{jet};#gamma_{e}", 3, -.5, 2.5, 3, -.5, 2.5, 3, -.5, 2.5 )
{

	for( int i = 0; i<nFiles; ++i )
		inputTree.Add( fileList[i] );
	event.setInput( inputTree );

	// Here the number of proceeded events will be stored. For plotting, simply use L*sigma/eventNumber
	eventNumbers.GetXaxis()->SetBinLabel(1,"Number of generated events");

	// Define one dimensional histograms
	hist1D["gMet"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["eMet"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["fMet"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["fMetUp"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["fMetDown"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["metFilters"] = TH1F("", ";met Filter number;Entries", 19, -.5, 18.5 );

	// Define two dimensional histograms
	hist2D["matchJet"] = TH2F("", "photon-jet matching;#DeltaR;p_{T, jet}/p_{T, #gamma}", 100, 0, 1, 100, 0, 4 );
	hist2D["matchJetFO"] = TH2F("", "photon-jet matching;#DeltaR;p_{T, jet}/p_{T, #gamma}", 100, 0, 1, 100, 0, 4 );
	hist2D["matchJetPt"] = TH2F("", "photon-jet matching;p_{T,#gamma};p_{T, jet}", 100, 0, 1000, 100, 0, 1000 );
	hist2D["matchJetPtFO"] = TH2F("", "photon-jet matching;p_{T,#gamma};p_{T, jet}", 100, 0, 1000, 100, 0, 1000 );
	hist2D["matchPhoton"] = TH2F("", ";#DeltaR;#Delta p_{T}/p_{T}", 100, 0, .5, 200, -2, 2 );
	hist2D["matchElectron"] = TH2F("", ";#DeltaR;#Delta p_{T}/p_{T}", 100, 0, .5, 200, -2, 2 );

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
		std::cerr << "ERROR: Could not read qcd weight file " << filename << std::endl;

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

	for( std::vector<const char*>::const_iterator it = triggerNames.begin();
			it != triggerNames.end(); ++it ) {
		for( susy::TriggerMap::const_iterator tm = event.hltMap.begin();
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
	int bin = qcdWeightHisto.FindBin( pt, ht_ );
	float error = qcdWeightHisto.GetBinError(bin);
	qcdWeight = qcdWeightHisto.GetBinContent(bin);
	qcdWeightUp = qcdWeight + error;
	qcdWeightDown = qcdWeight - error;
}

void TreeWriter::getPtFromMatchedJet( tree::Photon& myPhoton, bool isPhoton=false, bool isPhotonJet=false ) {
	/**
	 * \brief Takes jet p_T as photon p_T
	 *
	 * At first all jets with DeltaR < 0.3 (isolation cone) are searched.
	 * If several jets are found, take the one with the minimal pt difference
	 * compared to the photon. If no such jets are found, keep the photon_pt

	 * change _ptJet and matchedJetIndex
	 */
	myPhoton._ptJet = 0;
	myPhoton.matchedJetIndex = -1;
	std::vector<short> indices;

	for(std::vector<tree::Jet>::iterator jet = jets.begin();
			jet != jets.end(); ++jet) {

		float deltaR_ = myPhoton.DeltaR( *jet );
		float eRel = jet->pt / myPhoton.pt;

		// Fill in matching histograms
		if( isPhoton )
			hist2D["matchJet"].Fill( deltaR_, eRel, weight );
		if( isPhoton && deltaR_ < .2 )
			hist2D["matchJetPt"].Fill( myPhoton.pt, jet->pt, weight );
		if( isPhotonJet )
			hist2D["matchJetFO"].Fill( deltaR_, eRel, weight );
		if( isPhotonJet && deltaR_ < .2 )
			hist2D["matchJetPtFO"].Fill( myPhoton.pt, jet->pt, weight );

		// Define the selection criteria
		if( deltaR_ > 0.2  && eRel > 3 ) continue;
		jet->setMatch( tree::kJetAllPhoton );

		// If only one jet is found, we would be done here
		myPhoton.matchedJetIndex = std::distance( jets.begin(), jet );
		myPhoton._ptJet = jet->pt;

		// If more than one jet is found, we have to decide which jet to choose
		indices.push_back( myPhoton.matchedJetIndex );
	}// for jet

	// If more than one jet was found, we take the one nearer to the photon in pt
	if( indices.size() > 1 ) {
		float minPtDifferenz = 1E20; // should be very high

		for( std::vector<short>::const_iterator index = indices.begin();
				index != indices.end(); ++index ) {
			float ptDiff = std::abs(myPhoton.pt - jets.at(*index).pt);
			if ( ptDiff < minPtDifferenz ) {
				minPtDifferenz = ptDiff;
				myPhoton._ptJet = jets.at(*index).pt;
				myPhoton.matchedJetIndex = *index;
			}
		}
		jets.at( myPhoton.matchedJetIndex ).setMatch( tree::kJetPhoton );
	} else if( indices.size() == 1 ) {
		jets.at(0).setMatch( tree::kJetPhoton );
	} else if( loggingVerbosity > 1 )
		std::cout << "No matching jet found, do not change photon_pt." << std::endl;
}

void TreeWriter::fillJets() {
	/* Read the jets from susyEvent and save them to jet vector.
	 * All jets for HT calculation, and photon-jet matching are saved.
	 * This are not the final jets in the analysis.
	 * All jets are corrected
	 */
	jets.clear();
	tree::Jet jetToTree;

	std::vector<susy::PFJet> jetVector = event.pfJets.find("ak5chs")->second;
	for(std::vector<susy::PFJet>::const_iterator it = jetVector.begin();
			it != jetVector.end(); ++it) {

		TLorentzVector corrP4 = it->jecScaleFactors.at("L1FastL2L3") * it->momentum;

		if( std::abs(corrP4.Eta()) > 3 ) continue;
		if( corrP4.Pt() < 30 ) continue;
		if( !passLooseJetId( *it ) ) continue;

		jetToTree.matchInformation = 0;
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
	/* HT is sum jet + sum photon + sum photonJet + sum photonElectron
	 * For the jet sum, jets have to have a good Id or was used as match for a
	 * photon/photonJet/photonElectron.
	 * For the sum of photonObjects, the pt is only added in case the pt of the
	 * matched jet was not added.
	 */
	if( !splitting )
		return 0;

	float returnedHt = 0;
	for(std::vector<tree::Jet>::const_iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {

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

unsigned int TreeWriter::countGoodJets( bool clean ) {
	/* Count the number of good jets.
	 * They
	 * * have different pt and eta criteria as the jet collection
	 * * jetID
	 * * cleared of electorns, muons, photonObjects
	 */
	unsigned int number = 0;
	for(std::vector<tree::Jet>::iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {
		if( jet->pt < 30 || jet->eta > 2.5 ) continue;

		if( isAdjacentToParticles<tree::Particle>( *jet, electrons ) ) continue;
		if( isAdjacentToParticles<tree::Particle>( *jet, muons ) ) continue;
		if( clean ) {
			if( isAdjacentToParticles<tree::Photon>( *jet, photons ) ) continue;
			if( isAdjacentToParticles<tree::Photon>( *jet, photonElectrons ) ) continue;
			if( isAdjacentToParticles<tree::Photon>( *jet, photonJets ) ) continue;
		}
		jet->setMatch( tree::kJetCount );
		++number;
	}
	return number;
}
void TreeWriter::SetBranches( TTree& tree ) {
	/* For each tree, the branches have to be set
	 */
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
		if( printCascade ) {
			std::cout << "=================================================" << std::endl;
			for( susy::ParticleCollection::const_iterator it = event.genParticles.begin(); printCascade && it != event.genParticles.end(); ++it ){
				if( it->motherIndex == -1 ){
					printChildren( std::distance<susy::ParticleCollection::const_iterator>(event.genParticles.begin(), it ), event.genParticles );
				}
			}
			continue;
		}

		if ( event.isRealData )
			if ( !isGoodLumi() || !passTrigger()) continue;

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

		// For data, the weight is 1. Else take the pileup weight.
		weight = event.isRealData ? 1. : getPileUpWeight();

		// Particles needed for jet matching
		std::vector<tree::Particle> genQuarkLike;
		genQuarkLike.clear();

		// genParticles
		tree::Particle thisGenParticle;
		for( std::vector<susy::Particle>::const_iterator it = event.genParticles.begin(); it != event.genParticles.end(); ++it ) {

			// status 3: particles in matrix element
			// status 2: intermediate particles
			// status 1: final particles (but can decay in geant, etc)
			if( it->momentum.Pt() < 40 || it->status != 1) continue;

			thisGenParticle.pt = it->momentum.Pt();
			thisGenParticle.eta = it->momentum.Eta();
			thisGenParticle.phi = it->momentum.Phi();
			int pdgId = std::abs(it->pdgId);
			switch( pdgId ) {
				case 22: // photon
					genPhotons.push_back( thisGenParticle );
					break;
				case 11: // electron
					genElectrons.push_back( thisGenParticle );
					break;
			}
			if( pdgId < 7 // quarks
				|| (pdgId > 99 && pdgId < 6000) ) {// hadrons
				genQuarkLike.push_back( thisGenParticle );
			}
		}

		// electrons
		std::vector<susy::Electron> eVector = event.electrons["gsfElectrons"];
		for(std::vector<susy::Electron>::const_iterator it = eVector.begin(); it < eVector.end(); ++it) {
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
		for( std::vector<susy::Muon>::const_iterator it = mVector.begin(); it != mVector.end(); ++it) {
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

		// The jets have to be filled before looping over the photons and searching
		// for jet photon matches.
		fillJets();
		if( loggingVerbosity > 1 )
			std::cout << "Found " << jets.size() << " uncleaned jets" << std::endl;

		// photons
		std::vector<susy::Photon> photonVector = event.photons["photons"];
		for(std::vector<susy::Photon>::iterator it = photonVector.begin();
				it != photonVector.end(); ++it ) {

			if( std::abs( it->momentum.Eta() ) >= susy::etaGapBegin ) continue;

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

			//photon definition barrel
			bool isPhotonOrElectron = photonToTree.hadTowOverEm < 0.05
				&& photonToTree.sigmaIetaIeta < 0.012
				&& photonToTree.chargedIso < 2.6
				&& photonToTree.neutralIso < 3.5+0.04*photonToTree.pt
				&& photonToTree.photonIso < 1.3+0.005*photonToTree.pt;

			bool isPhoton = isPhotonOrElectron && !photonToTree.pixelseed;
			bool isPhotonElectron = isPhotonOrElectron && photonToTree.pixelseed;

			// photonJet definition
			bool isPhotonJet = !photonToTree.pixelseed
				&& photonToTree.hadTowOverEm < 0.05
				&& photonToTree.sigmaIetaIeta < 0.012
				&& photonToTree.chargedIso < 26 && photonToTree.chargedIso > 0.26
				&& photonToTree.neutralIso < 35+0.4*photonToTree.pt && photonToTree.neutralIso > 0.35+0.004*photonToTree.neutralIso
				&& photonToTree.photonIso < 13+0.05*photonToTree.pt && photonToTree.photonIso > 0.13+0.0005*photonToTree.pt;

			// print photon information
			if( loggingVerbosity > 2 ) {
				if( isPhoton )         std::cout << " photon pT = " << photonToTree.pt << std::endl;
				if( isPhotonElectron ) std::cout << " photonElectron pT = " << photonToTree.pt << std::endl;
				if( isPhotonJet )      std::cout << " photonJet pT = " << photonToTree.pt << std::endl;
			}

			// Fill matching histograms only for photon-like objects
			if( splitting && !isPhotonOrElectron && !isPhotonJet ) continue;

			if( matchLorentzToGenVector( it->momentum, genPhotons, &hist2D["matchPhoton"], .1 ) )
				photonToTree.setGen( tree::kGenPhoton );
			if( matchLorentzToGenVector( it->momentum, genElectrons, &hist2D["matchElectron"], .1 ) )
				photonToTree.setGen( tree::kGenElectron );
			if( matchLorentzToGenVector( it->momentum, genQuarkLike, NULL, .3 ) )
				photonToTree.setGen( tree::kGenJet );

			if( matchLorentzToGenVector( it->momentum, electrons, NULL, .3 ) ||
					matchLorentzToGenVector( it->momentum, muons, NULL, .3 ) )
				photonToTree.setGen( tree::kNearLepton );

			getPtFromMatchedJet( photonToTree, isPhoton, isPhotonJet );
			if( loggingVerbosity > 2 )
				std::cout << "  ->jet pT = " << photonToTree._ptJet << std::endl;

			// If no jet is found for a loose photon, the photon is rejected
			if( isPhotonJet && !photonToTree._ptJet ) continue;

			if( photonToTree.ptJet() < photonPtThreshold ) continue;

			if( splitting ) {
					if( isPhotonOrElectron ) {
						if( photonToTree.pixelseed )
							photonElectrons.push_back( photonToTree );
						else
							photons.push_back( photonToTree );
					} else if ( isPhotonJet )
						photonJets.push_back( photonToTree );
			} else // no splitting, put everything in the vector 'photons'
				photons.push_back( photonToTree );
		}
		std::sort( photons.begin(), photons.end(), tree::EtGreater );
		std::sort( photonElectrons.begin(), photonElectrons.end(), tree::EtGreater );
		std::sort( photonJets.begin(), photonJets.end(), tree::EtGreater );
		if( loggingVerbosity > 1 )
			std::cout << "Found " << photons.size() << " photons, "
					<< photonJets.size() << " photon_{jets} and "
					<< photonElectrons.size() << " photon electrons." << std::endl;
		nPhotons.Fill( photons.size(), photonJets.size(), photonElectrons.size() );

		// filter out events with no photons
		if( !photons.size() && !photonJets.size() && !photonElectrons.size() ) continue;

		// met
		met = event.metMap["pfMet"].met();
		type0met = event.metMap["pfType01CorrectedMet"].met();
		type1met = event.metMap["pfType1CorrectedMet"].met();
		if( loggingVerbosity > 2 )
			std::cout << " met = " << met << std::endl;

		ht = getHt();
		nGoodJets = countGoodJets( splitting );

		//if( event.passMetFilters() ) continue;
		if( splitting && hadronicSelection && ( nGoodJets < 2 || ht < 500 ) ) continue;

		fillMetFilterBitHistogram( hist1D.at("metFilters"), event.metFilterBit );
		if( !event.passMetFilters() ) continue;

		if( splitting ) {
			// Assing event to the leading photonObject
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
		for( std::map<std::string, TH2F>::const_iterator it = hist2D.begin();
				it!= hist2D.end(); ++it )
			it->second.Write();
	}

	// If running over signal scans, the mass point information is appended to
	// the histogram name.
	TPRegexp expFilename( ".*/tree_([0-9]+_[0-9]+)_375.root" ); // eg. /path/to/mc/tree_1200_1220_375.root
	TObjArray *arr = expFilename.MatchS( inputTree.GetCurrentFile()->GetName() );

	std::string histoNameAppendix = "";
	if( arr->GetLast() >0 )
		histoNameAppendix = (std::string)(((TObjString *)arr->At(1))->GetString());
	else if( loggingVerbosity > 0 )
		std::cout << "Could not extract grid parameters from filename." << std::endl;

	// Append the signal information to the histogram name
	for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
			it!= hist1D.end(); ++it ) {
		it->second.SetName( (it->second.GetName() + histoNameAppendix ).c_str() );
		it->second.Write();
	}
}

