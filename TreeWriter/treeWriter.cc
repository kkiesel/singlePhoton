#include "treeWriter.h"
#include "printCascade.h"

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

tree::electronWorkingPoints getElectronWorkingPoint ( const susy::Electron& electron, const susy::Event& event ) {
	/** Definition of all electron working points.
	 *
	 * See https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
	 * for more information.
	 */
	if( electron.momentum.Pt() > 1e6 )
		return tree::kNoElectron; // spike rejection

	// list and compute all variables
	float fabsdEtaIn = fabs(electron.deltaEtaSuperClusterTrackAtVtx);
	float fabsdPhiIn = fabs(electron.deltaPhiSuperClusterTrackAtVtx);
	// electron.sigmaIetaIeta
	// electron.hcalOverEcalBc
	float d0 = fabs( electron.gsfTrack->d0( electron.vertex ) );
	float dZ = fabs( electron.gsfTrack->dz( electron.vertex ) );
	float fabsInvDiff = fabs( 1./electron.ecalEnergy - 1./electron.trackMomentumAtVtx.Pt() );
	float iso = ( electron.chargedHadronIso +
		std::max(electron.neutralHadronIso+electron.photonIso -
		effectiveAreaElectron(electron.momentum.Eta())*event.rho25, (float)0. ))
		/ electron.momentum.Pt();
	// electron.passConversionVeto
	// electron.nMissingHits
	float eta = std::abs(electron.superCluster->position.Eta());
	bool isBarrel = eta <= 1.479;
	bool isEndcap = eta > 1.479 && eta < 2.5;

	if( ( isBarrel
			&& fabsdEtaIn < 0.004
			&& fabsdPhiIn < 0.03
			&& electron.sigmaIetaIeta < 0.01
			&& electron.hcalOverEcalBc < 0.12
			&& d0 < 0.02
			&& dZ < 0.1
			&& fabsInvDiff < 0.05
			&& iso < 0.1
			&& electron.passConversionVeto
			&& electron.nMissingHits == 0
		) || ( isEndcap
			&& fabsdEtaIn < 0.005
			&& fabsdPhiIn < 0.02
			&& electron.sigmaIetaIeta < 0.03
			&& electron.hcalOverEcalBc < 0.1
			&& d0 < 0.02
			&& dZ < 0.1
			&& fabsInvDiff < 0.05
			&& iso < 0.1
			&& ( iso < 0.07 || electron.momentum.Pt() > 20 )
			&& electron.passConversionVeto
			&& electron.nMissingHits == 0
		) )
		return tree::kTightElectron;

	if( ( isBarrel
			&& fabsdEtaIn < 0.004
			&& fabsdPhiIn < 0.06
			&& electron.sigmaIetaIeta < 0.01
			&& electron.hcalOverEcalBc < 0.12
			&& d0 < 0.02
			&& dZ < 0.1
			&& fabsInvDiff < 0.05
			&& iso < 0.15
			&& electron.passConversionVeto
			&& electron.nMissingHits <= 1
		) || ( isEndcap
			&& fabsdEtaIn < 0.007
			&& fabsdPhiIn < 0.03
			&& electron.sigmaIetaIeta < 0.03
			&& electron.hcalOverEcalBc < 0.1
			&& d0 < 0.02
			&& dZ < 0.1
			&& fabsInvDiff < 0.05
			&& iso < 0.15
			&& ( iso < 0.10 || electron.momentum.Pt() > 20 )
			&& electron.passConversionVeto
			&& electron.nMissingHits <= 1
		) )
		return tree::kMediumElectron;

	if( ( isBarrel
			&& fabsdEtaIn < 0.007
			&& fabsdPhiIn < 0.15
			&& electron.sigmaIetaIeta < 0.01
			&& electron.hcalOverEcalBc < 0.12
			&& d0 < 0.02
			&& dZ < 0.2
			&& fabsInvDiff < 0.05
			&& iso < 0.15
			&& electron.passConversionVeto
			&& electron.nMissingHits <= 1
		) || ( isEndcap
			&& fabsdEtaIn < 0.009
			&& fabsdPhiIn < 0.1
			&& electron.sigmaIetaIeta < 0.03
			&& electron.hcalOverEcalBc < 0.1
			&& d0 < 0.02
			&& dZ < 0.2
			&& fabsInvDiff < 0.05
			&& iso < 0.15
			&& ( iso < 0.10 || electron.momentum.Pt() > 20 )
			&& electron.passConversionVeto
			&& electron.nMissingHits <= 1
		) )
		return tree::kLooseElectron;

	if( ( isBarrel
			&& fabsdEtaIn < 0.007
			&& fabsdPhiIn < 0.8
			&& electron.sigmaIetaIeta < 0.01
			&& electron.hcalOverEcalBc < 0.15
			&& d0 < 0.04
			&& dZ < 0.2
			&& iso < 0.15
		) || ( isEndcap
			&& fabsdEtaIn < 0.01
			&& fabsdPhiIn < 0.7
			&& electron.sigmaIetaIeta < 0.03
			&& d0 < 0.04
			&& dZ < 0.2
			&& iso < 0.15
		) )
		return tree::kVetoElectron;

	return tree::kNoElectron;
}

bool isLooseJet( const susy::PFJet& jet ) {
	/**
	 * \brief Apply loose cut on jets.
	 *
	 * See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_7_TeV_data_a
	 * for more information.
	 */
	double energy = jet.momentum.E();
	//double energy = jet.chargedHadronEnergy+jet.neutralHadronEnergy + jet.photonEnergy+jet.electronEnergy+jet.muonEnergy+jet.HFHadronEnergy+jet.HFEMEnergy;
	bool debug = false;
	// debug = true;
	if( debug )
		std::cout << "Jet Id with pt = " << jet.jecScaleFactors.at("L1FastL2L3") * jet.momentum.Pt()
			<< "\nNeutralHadron = " << (jet.neutralHadronEnergy+jet.HFHadronEnergy) / energy
			<< "\nNeutralEM = " << jet.neutralEmEnergy / energy
			<< "\nnConst = " << (int)jet.nConstituents
			<< "\neta = " << std::abs(jet.momentum.Eta())
			<< "\nchargedHadron = " << jet.chargedHadronEnergy / energy
			<< "\nnCharged = " << (int)jet.chargedMultiplicity
			<< "\nchargedEm = " << jet.chargedEmEnergy / energy
			<< std::endl;

	return (jet.neutralHadronEnergy+jet.HFHadronEnergy) / energy < 0.99
			&& jet.neutralEmEnergy / energy < 0.99
			&& jet.nConstituents > 1
			&& ( std::abs(jet.momentum.Eta()) >= 2.4
				|| ( jet.chargedHadronEnergy / energy > 0
					&& jet.chargedMultiplicity > 0
					&& jet.chargedEmEnergy / energy < 0.99 ) );
}

bool goodVertex( const susy::Vertex& vtx ) {
	/** Definition of a good vertex. Returns true if the vertex is good.
	 */
	return (!vtx.isFake() &&
		vtx.ndof >= 4 &&
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

unsigned int nTrackPrimaryVertex( const std::vector<susy::Vertex>& vertexVector ) {
	/* Tracks coming from the first good vertex */
	for( std::vector<susy::Vertex>::const_iterator vtx = vertexVector.begin();
			vtx != vertexVector.end(); ++vtx ) {
		if( goodVertex( *vtx ) )
			return vtx->tracksSize;
	}
	// if no valid vertex was found
	return 0;
}

template <typename VectorClass>
int indexOfnearestParticle( const tree::Particle& thisParticle, const std::vector<VectorClass>& particleVector,
		float deltaR_=.3, float ptRelMin_=-1e6, float ptRelMax_=1e6, TH2F* hist=NULL ) {
	/* Compare 'thisParticle' to each particle in the vector
	 *
	 * If a particle in the vector satisfies the cut values, its index is returned.
	 * If no particle if found, -1 is returned.
	 * If serveral particles satisfy the requirements, the nearest in deltaR is
	 * returned.
	 */

	int index = -1; //default value
	std::map<float,int> map_dr_i;
	float dR, relPt;
	for( typename std::vector<VectorClass>::const_iterator it = particleVector.begin();
			it != particleVector.end(); ++it ) {

		relPt = it->pt / thisParticle.pt;
		dR = thisParticle.DeltaR( *it );

		if( hist ) hist->Fill( dR, relPt );

		if( dR > deltaR_ ||  relPt > ptRelMax_ || relPt < ptRelMin_ ) continue;
		index = std::distance( particleVector.begin(), it );
		map_dr_i[dR] = index;
	}

	// search nearest if several particles found
	if( map_dr_i.size() > 1 ) {
		index = map_dr_i.begin()->second;
	}

	return index;
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

void tryFill( std::map< std::string, TH1F >& histMap, const std::string& histname, const std::string& appendix, float var, float weight=1 ) {
	if ( histMap.find( histname+appendix ) == histMap.end() ) {
		histMap[ histname+appendix ] = *((TH1F*) histMap[ histname ].Clone( (histname+appendix).c_str() ));
		histMap[ histname+appendix ].Reset( "ICESM" );
	}
	histMap[ histname+appendix ].Fill( var, weight );
}


///////////////////////////////////////////////////////////////////////////////
// Here the class implementation begins ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

TreeWriter::TreeWriter( int nFiles, char** fileList, std::string const& outputName) :
	reportEvery(200000),
	processNEvents(-1),
	loggingVerbosity(0),
	inputTree("susyTree"),
	event(),
	outFile( outputName.c_str(), "recreate" ),
	photonTree("photonTree","Tree for single photon analysis"),
	photonElectronTree("photonElectronTree","Tree for single photon analysis"),
	photonJetTree("photonJetTree","Tree for single photon analysis"),
	eventNumbers("eventNumbers", "Histogram containing number of generated events", 1, 0, 1),
	nPhotons("nPhotons", ";#gamma_{tight};#gamma_{loose};#gamma_{pixel}", 3, -.5, 2.5, 3, -.5, 2.5, 3, -.5, 2.5 )
{

	for( int i = 0; i<nFiles; ++i )
		inputTree.Add( fileList[i] );
	event.setInput( inputTree );

	// Here the number of proceeded events will be stored. For plotting, simply use L*sigma/eventNumber
	eventNumbers.GetXaxis()->SetBinLabel(1,"Number of generated events");

	// Define one dimensional histograms
	hist1D["gMet"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["nGen"] = TH1F("", ";met;", 1, 0, 1 );
	hist1D["gMetPuUp"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["gMetPuDown"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["gMetJesUp"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["gMetJesDown"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["eMet"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["fMet"] = TH1F("", ";met;", 60, 0, 600 );
	hist1D["fMetError"] = TH1F("", ";met;", 60, 0, 600 );

	hist1D["metFilters"] = TH1F("", ";met Filter number;Entries", susy::nMetFilters+1, .5, susy::nMetFilters+1.5 );
	hist1D["metFilters"].Fill( 0., 0. ); // Allows the histograms to be merged

	hist1D["gHt"] = TH1F("", ";H_{T} [GeV];Entries", 200, 0, 2000 );
	hist1D["gNJets"] = TH1F("", ";n_{Jets};Entries", 10, -.5, 9.5 );
	hist1D["gPt"] = TH1F("", ";p_{T^{*}};Entries", 200, 0, 2000 );

	// Define two dimensional histograms
	hist2D["matchPhotonToJet"]         = TH2F("", "photon-jet matching;#DeltaR;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 1, 100, 0, 4 );
	hist2D["matchPhotonJetToJet"]      = TH2F("", "photon-jet matching;#DeltaR;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 1, 100, 0, 4 );
	hist2D["matchPhotonElectronToJet"] = TH2F("", "photon-jet matching;#DeltaR;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 1, 100, 0, 4 );

	hist2D["matchGenPhoton"]   = TH2F("", ";#DeltaR;p_{T}^{gen} / p_{T}", 1000, 0, .5, 200, 0, 2 );
	hist2D["matchGenElectron"] = TH2F("", ";#DeltaR;p_{T}^{gen} / p_{T}", 1000, 0, .5, 200, 0, 2 );

	std::string histoNameAppendix = "";
	if( runType == kGMSB ) {
		// If running over signal scans, the mass point information is appended to
		// the histogram name.
		TPRegexp expFilename( ".*/tree_([0-9]+_[0-9]+)_375.root" ); // eg. /path/to/mc/tree_1200_1220_375.root
		TObjArray *arr = expFilename.MatchS( inputTree.GetCurrentFile()->GetName() );

		if( arr->GetLast() >0 )
			histoNameAppendix = (std::string)(((TObjString *)arr->At(1))->GetString());
		else if( loggingVerbosity > 0 )
			std::cout << "Could not extract grid parameters from filename." << std::endl;
	}

	// Set the keyName as histogram name for one and two dimensional histograms
	for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
			it!= hist1D.end(); ++it ) {
		it->second.SetName( (it->second.GetName() + it->first + histoNameAppendix ).c_str() );
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

eventType TreeWriter::whichEventType( std::vector<tree::Photon>& photons_,
			std::vector<tree::Photon>& photonElectrons_,
			std::vector<tree::Photon>& photonJets_ ) {

	float gPt = photons_.size()         ? photons_.at(0).pt         : 0;
	float ePt = photonElectrons_.size() ? photonElectrons_.at(0).pt : 0;
	float fPt = photonJets_.size()      ? photonJets_.at(0).pt      : 0;

	if( gPt && gPt > fPt && gPt > ePt )
		return kPhotonEvent;
	else if( ePt && ePt > fPt && ePt > gPt )
		return kElectronEvent;
	else if( fPt && fPt > ePt && fPt > gPt )
		return kJetEvent;
	else{
		std::cerr << "This should not happen" << std::endl;
		return kPhotonEvent;
	}
}


float TreeWriter::getPileUpWeight(){
	/**
	 * If a pileup weight histogram has been added, the pile-up weight for the
	 * current event is computed.
	 */

	float thisWeight = 0;
	for( susy::PUSummaryInfoCollection::const_iterator iBX = event.pu.begin();
			iBX != event.pu.end(); ++iBX) {
		if (iBX->BX == 0) { // find bunch crossing for this event
			float trueNumInteractions = iBX->trueNumInteractions;
			thisWeight = pileupHisto.GetBinContent( pileupHisto.FindBin( trueNumInteractions ) );
			weightPuUp = pileupHistoUp.GetBinContent( pileupHistoUp.FindBin( trueNumInteractions ) );
			weightPuDown = pileupHistoDown.GetBinContent( pileupHistoDown.FindBin( trueNumInteractions ) );
			break;
		}
	}

	if( loggingVerbosity > 2 )
		std::cout << "Pile-up weight = " << thisWeight << std::endl;
	return thisWeight;
}

void TreeWriter::getQcdWeights( float pt, float ht_, float & qcdWeight, float & qcdWeightError ){
	int bin = qcdWeightHisto.FindBin( pt, ht_ );
	qcdWeight = qcdWeightHisto.GetBinContent(bin);
	qcdWeightError = qcdWeightHisto.GetBinError(bin);
}

void TreeWriter::fillGenParticles() {
	genPhotons.clear();
	genElectrons.clear();

	// genParticles
	tree::Particle thisGenParticle;
	for( susy::ParticleCollection::const_iterator it = event.genParticles.begin();
			it != event.genParticles.end(); ++it ) {

		// status 3: particles in matrix element
		// status 2: intermediate particles
		// status 1: final particles (but can decay in geant, etc)
		if( it->momentum.Pt() < 40 || it->status != 1) continue;

		thisGenParticle.pt = it->momentum.Pt();
		thisGenParticle.eta = it->momentum.Eta();
		thisGenParticle.phi = it->momentum.Phi();
		thisGenParticle.bitFlag = 0;
		int pdgId = std::abs(it->pdgId);
		switch( pdgId ) {
			case 22: // photon
				genPhotons.push_back( thisGenParticle );
				break;
			case 11: // electron
				genElectrons.push_back( thisGenParticle );
				break;
		}
	}

	if( loggingVerbosity > 1 )
		std::cout << "Found " << genPhotons.size() << " generated photons and "
		<< genElectrons.size() << " generated electrons" << std::endl;
}

void TreeWriter::fillLeptons() {
	electrons.clear();
	muons.clear();
	tree::Particle leptonToTree;

	// electrons
	std::vector<susy::Electron> eVector = event.electrons["gsfElectrons"];
	for(std::vector<susy::Electron>::const_iterator it = eVector.begin(); it < eVector.end(); ++it) {
		if( it->momentum.Pt() < 15 || std::abs(it->momentum.Eta()) > 2.6 )
			continue;
		tree::electronWorkingPoints wp = getElectronWorkingPoint( *it, event );
		if( wp == tree::kNoElectron )
			continue;
		leptonToTree.pt = it->momentum.Pt();
		leptonToTree.eta = it->momentum.Eta();
		leptonToTree.phi = it->momentum.Phi();
		leptonToTree.bitFlag = 0;
		leptonToTree.setStatus( wp );
		if( indexOfnearestParticle<tree::Photon>( leptonToTree, photons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Photon>( leptonToTree, photonElectrons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Photon>( leptonToTree, photonJets, .3 ) > -1 ) continue;
		electrons.push_back( leptonToTree );
	}
	if( loggingVerbosity > 1 )
		std::cout << "Found " << electrons.size() << " electrons" << std::endl;

	// muons
	std::vector<susy::Muon> mVector = event.muons["muons"];
	for( std::vector<susy::Muon>::const_iterator it = mVector.begin(); it != mVector.end(); ++it) {
		// see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
		if( it->momentum.Pt() < 15 || std::abs(it->momentum.Eta()) > 2.6 || !(it->isPFMuon() && (it->isGlobalMuon() || it->isTrackerMuon())) )
			continue;
		leptonToTree.pt = it->momentum.Et();
		leptonToTree.eta = it->momentum.Eta();
		leptonToTree.phi = it->momentum.Phi();
		leptonToTree.bitFlag = 0;
		if( indexOfnearestParticle<tree::Photon>( leptonToTree, photons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Photon>( leptonToTree, photonElectrons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Photon>( leptonToTree, photonJets, .3 ) > -1 ) continue;
		muons.push_back( leptonToTree );
	}
	if( loggingVerbosity > 1 )
		std::cout << "Found " << muons.size() << " muons" << std::endl;
}

void TreeWriter::fillJets( int jecScale=0 ) {
	/* Read the jets from susyEvent and save them to jet vector.
	 * All jets for HT calculation, and photon-jet matching are saved.
	 * This are not the final jets in the analysis.
	 * All jets are corrected
	 */
	jets.clear();
	tree::Jet jetToTree;

	/*std::string jecDir("../../CMSSW/CMSSW_5_3_8_patch3/src/SUSYPhotonAnalysis/SusyNtuplizer/jec/");
	FactorizedJetCorrector jecCHS("L1FastJet:L2Relative:L3Absolute:L2Relative",
		jecDir + "FT_53_V21_AN5_L1FastJet_AK5PFchs.txt:" +
		jecDir + "FT_53_V21_AN5_L2Relative_AK5PFchs.txt:" +
		jecDir + "FT_53_V21_AN5_L3Absolute_AK5PFchs.txt:" +
		jecDir + "FT_53_V21_AN5_L2RelativeL3AbsoluteResidual_AK5PFchs.txt" );
	*/

	std::vector<susy::PFJet> jetVector = event.pfJets.find("ak5chs")->second;
	for(std::vector<susy::PFJet>::const_iterator it = jetVector.begin();
			it != jetVector.end(); ++it) {

		TLorentzVector corrP4 = it->jecScaleFactors.at("L1FastL2L3") * it->momentum;
		corrP4 *= (1 + jecScale*it->jecUncertainty );

		if( std::abs(corrP4.Eta()) > 3 ) continue;
		if( corrP4.Pt() < 30 ) continue;

		/*std::cout << "\nNew Jet with " << it->momentum.Pt()<<"\n";
		std::cout << "area = " << it->jetArea << " rho = " << event.rho << std::endl;
		std::cout << "nTuple scale = " << it->jecScaleFactors.at("L1FastL2L3") <<std::endl;
		jecCHS.setJetEta(it->momentum.Eta());
		jecCHS.setJetPt(it->momentum.Pt());
		jecCHS.setJetA(it->jetArea);
		jecCHS.setRho(event.rho);
		std::vector<float> subcorr = jecCHS.getSubCorrections();
		std::cout << "uncorrected pt = " << it->momentum.Pt() << "   1" << std::endl;
		std::cout << "L1corrected pt = " << it->momentum.Pt()*subcorr.at(0) << "   " << subcorr.at(0) <<  std::endl;
		std::cout << "L2corrected pt = " << it->momentum.Pt()*subcorr.at(1) << "   " << subcorr.at(1) <<  std::endl;
		std::cout << "L3corrected pt = " << it->momentum.Pt()*subcorr.at(2) << "   " << subcorr.at(2) <<  std::endl;
		std::cout << "LRcorrected pt = " << it->momentum.Pt()*subcorr.at(3) << "   " << subcorr.at(3) <<  std::endl;
		*/

		jetToTree.bitFlag = 0;
		if( isLooseJet( *it ) )
			jetToTree.setStatus( tree::kJetId );


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

	if( loggingVerbosity > 1 )
			std::cout << "Found " << jets.size() << " uncleaned jets" << std::endl;
}

TVector3 TreeWriter::getRecoilVector( eventType eType ) const {
	TVector3 sum,adding;

	for(std::vector<tree::Jet>::const_iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {

		if( jet->pt < 30 || std::abs(jet->eta) > 3 ) continue;
		if( eType == kPhotonEvent && jet->DeltaR( photons[0])<0.5 ) continue;
		if( eType == kJetEvent && jet->DeltaR( photonJets[0])<0.5 ) continue;
		if( eType == kElectronEvent && jet->DeltaR( photonElectrons[0])<0.5 ) continue;

		adding.SetPtEtaPhi( jet->pt, jet->eta, jet->phi );
		sum += adding;
	}

	return sum;
}

TVector3 TreeWriter::getMhtVector() const {
	TVector3 sum,adding;

	for(std::vector<tree::Jet>::const_iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {

		if( jet->pt < 40 || std::abs(jet->eta) > 3. ) continue;
		adding.SetPtEtaPhi( jet->pt, jet->eta, jet->phi );
		sum += adding;
	}

	for( std::vector<tree::Photon>::const_iterator photon = photons.begin();
			photon != photons.end(); ++photon ) {
		if( photon->_ptJet == 0 ) {
			adding.SetPtEtaPhi( photon->pt, photon->eta, photon->phi );
			sum += adding;
		}
	}
	for( std::vector<tree::Photon>::const_iterator photon = photonJets.begin();
			photon != photonJets.end(); ++photon ) {
		if( photon->_ptJet == 0 ) {
			adding.SetPtEtaPhi( photon->pt, photon->eta, photon->phi );
			sum += adding;
		}
	}
	for( std::vector<tree::Photon>::const_iterator photon = photonElectrons.begin();
			photon != photonElectrons.end(); ++photon ) {
		if( photon->_ptJet == 0 ){
			adding.SetPtEtaPhi( photon->pt, photon->eta, photon->phi );
			sum += adding;
		}
	}
	return sum;
}

float TreeWriter::getHt() const {
	/* HT is sum jet + sum photon + sum photonJet + sum photonElectron
	 * For the jet sum, jets have to have a good Id or was used as match for a
	 * photon/photonJet/photonElectron.
	 * For the sum of photonObjects, the pt is only added in case the pt of the
	 * matched jet was not added.
	 */

	float returnedHt = 0;
	for(std::vector<tree::Jet>::const_iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {

		if( jet->pt < 40 || std::abs(jet->eta) > 3. ) continue;
		if( !jet->isStatus( tree::kJetId ) ) continue;
		//std::cout << " add jet to HT " << jet->pt << std::endl;

		returnedHt += jet->pt;
	}

	for( std::vector<tree::Photon>::const_iterator photon = photons.begin();
			photon != photons.end(); ++photon ) {
		if( photon->_ptJet == 0 ) {
			//std::cout << " add photon to HT " << photon->pt << std::endl;
			returnedHt += photon->pt;
		}
	}
	for( std::vector<tree::Photon>::const_iterator photon = photonJets.begin();
			photon != photonJets.end(); ++photon ) {
		if( photon->_ptJet == 0 ) {
			//std::cout << " add photonJet to HT " << photon->pt << std::endl;
			returnedHt += photon->pt;
		}
	}
	for( std::vector<tree::Photon>::const_iterator photon = photonElectrons.begin();
			photon != photonElectrons.end(); ++photon ) {
		if( photon->_ptJet == 0 ) {
			//std::cout << " add photonElectron to HT " << photon->pt << std::endl;
			returnedHt += photon->pt;
		}
	}
	return returnedHt;
}

unsigned int TreeWriter::countGoodJets() {
	/* Count the number of good jets.
	 * They
	 * * have different pt and eta criteria as the jet collection
	 * * jetID
	 * * cleared of electorns, muons, photonObjects
	 */
	unsigned int number = 0;
	for(std::vector<tree::Jet>::iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {
		if( !jet->isStatus( tree::kJetId )) continue;
		if( jet->pt < 30 || std::abs(jet->eta) > 2.5 ) continue;

		if( indexOfnearestParticle<tree::Particle>( *jet, electrons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Particle>( *jet, muons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Photon>( *jet, photons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Photon>( *jet, photonElectrons, .3 ) > -1 ) continue;
		if( indexOfnearestParticle<tree::Photon>( *jet, photonJets, .3 ) > -1 ) continue;
		jet->setStatus( tree::kJetCount );
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
	tree.Branch("metSig", &metSig, "metSig/F");
	tree.Branch("met", &met, "met/F");
	tree.Branch("metPhi", &metPhi, "metPhi/F");
	tree.Branch("metShiftxy", &metShiftxy, "metShiftxy/F");
	tree.Branch("metShiftxyPhi", &metShiftxyPhi, "metShiftxyPhi/F");
	tree.Branch("met01corr", &met01corr, "met01corr/F");
	tree.Branch("met01corrPhi", &met01corrPhi, "met01corrPhi/F");
	tree.Branch("recoil", &recoil, "recoil/F");
	tree.Branch("recoilPhi", &recoilPhi, "recoilPhi/F");
	tree.Branch("mht", &mht, "mht/F");
	tree.Branch("mhtPhi", &mhtPhi, "mhtPhi/F");
	tree.Branch("ht", &ht, "ht/F");
	tree.Branch("weight", &weight, "weight/F");
	tree.Branch("nVertex", &nVertex, "nVertex/I");
	tree.Branch("nGoodJets", &nGoodJets, "nGoodJets/i");
	tree.Branch("nTracksPV", &nTracksPV, "nTracksPV/i");
	tree.Branch("runNumber", &runNumber, "runNumber/i");
	tree.Branch("eventNumber", &eventNumber, "eventNumber/i");
	tree.Branch("luminosityBlockNumber", &luminosityBlockNumber, "luminosityBlockNumber/i");
}

void TreeWriter::Loop( int jetScale ) {
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

	// clear all 1d histograms
	for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
			it!= hist1D.end(); ++it ) {
		it->second.Reset("ICES M");
	}

	for (long jentry=0; jentry < processNEvents; ++jentry) {
		event.getEntry(jentry);

		// Just for testing purpose (leave this uncommented)
		//if( event.eventNumber != 7302527 ) continue; loggingVerbosity = 5;
		if ( loggingVerbosity>1 || jentry%reportEvery==0 )
			std::cout << jentry << " / " << processNEvents << std::endl;

		// Uncomment this to just print the cascade on console
		//printCascade( event.genParticles ); continue;

		// for the simpified model, get signal point by looking at generated particles
		std::stringstream signalPointStringStream;
		if( runType == kSimplifiedModel ) {
			// get signal point from generated particles
			int mGluino=0, mLSP=0;
			for( susy::ParticleCollection::const_iterator it = event.genParticles.begin();
					it != event.genParticles.end(); ++it ) {
				if( it->status == 3 ) { // only particles from matrix element

					// neutralino mass rounded to 25, 75, 125, ..., starting at 25 GeV in steps of 50 GeV
					if( std::abs(it->pdgId) == 1000023 ) mLSP = -25+50*round((it->momentum.M()+25)/50);
					// gluino mass rounded to multiples of 50
					if( it->pdgId == 1000021 ) mGluino = 50*round(it->momentum.M()/50);

				}
			} // for generated particles
		signalPointStringStream << mGluino << "_" << mLSP;
		} // got signalString for simplified models
		std::string signalPointString = signalPointStringStream.str();
		tryFill( hist1D, "nGen", signalPointString, 0 );

		if ( event.isRealData )
			if ( !isGoodLumi() || !passTrigger()) continue;

		// vertices
		nVertex = numberOfGoodVertexInCollection( event.vertices );
		if( !nVertex ) continue;
		nTracksPV = nTrackPrimaryVertex( event.vertices );
		if( loggingVerbosity > 2 )
			std::cout << " nTracksPV = " << nTracksPV << std::endl;

		runNumber = event.runNumber;
		eventNumber = event.eventNumber;
		luminosityBlockNumber = event.luminosityBlockNumber;

		// For data, the weight is 1. Else take the pileup weight.
		weight = event.isRealData ? 1. : getPileUpWeight();

		// The jets have to be filled before looping over the photons and searching
		// for jet photon matches. The collections are not cross cleaned.
		fillJets( jetScale );
		fillGenParticles();

		// met
		susy::MET pfMet = event.metMap["pfMet"];
		met = pfMet.met();
		metSig = pfMet.significance;
		metPhi = pfMet.mEt.Phi();
		if( loggingVerbosity > 2 )
			std::cout << " met = " << met << std::endl;

		// subtract z boson
		susy::MET pfMetxy = event.metMap["pfMet"];
		for( susy::ParticleCollection::const_iterator it = event.genParticles.begin(); it != event.genParticles.end(); ++it ) {
			if( it->status != 3 || it->pdgId != 23) continue;
			TVector2 zBoson;
			zBoson.SetMagPhi( it->momentum.Pt(), it->momentum.Phi() );
			pfMetxy.mEt += zBoson;
		}
		metShiftxy = pfMetxy.met();
		metShiftxyPhi = pfMetxy.mEt.Phi();

		met01corr = event.metMap["pfType01CorrectedMet"].met();
		met01corrPhi = event.metMap["pfType01CorrectedMet"].mEt.Phi();

		photons.clear();
		photonJets.clear();
		photonElectrons.clear();
		// photons
		std::vector<susy::Photon> photonVector = event.photons["photons"];
		for(std::vector<susy::Photon>::iterator it = photonVector.begin();
				it != photonVector.end(); ++it ) {

			if( std::abs( it->momentum.Eta() ) > susy::etaGapBegin ) continue;

			photonToTree.chargedIso = chargedHadronIso_corrected(*it, event.rho);
			photonToTree.neutralIso = neutralHadronIso_corrected(*it, event.rho);
			photonToTree.photonIso = photonIso_corrected(*it, event.rho);
			photonToTree.pt = it->momentum.Pt();
			photonToTree.eta = it->momentum.Eta();
			photonToTree.phi = it->momentum.Phi();
			photonToTree.r9 = it->r9;
			photonToTree.sigmaIetaIeta = it->sigmaIetaIeta;
			photonToTree.sigmaIphiIphi = it->sigmaIphiIphi;
			photonToTree.hadTowOverEm = it->hadTowOverEm;
			photonToTree.pixelseed = it->nPixelSeeds;
			photonToTree.conversionSafeVeto = it->passelectronveto;
			photonToTree.bitFlag = 0;

			int jetIndex = indexOfnearestParticle<tree::Jet>( photonToTree, jets, .2, .8, 3 );
			photonToTree._ptJet = jetIndex>-1 ? jets.at(jetIndex).pt : 0;
			photonToTree._etaJet = jetIndex>-1 ? jets.at(jetIndex).eta : 0;
			photonToTree._phiJet = jetIndex>-1 ? jets.at(jetIndex).phi : 0;


			//photon definition barrel
			bool isPhotonOrElectron =
				( std::abs(photonToTree.eta) <= susy::etaGapBegin
					&& photonToTree.hadTowOverEm < 0.05
					&& photonToTree.sigmaIetaIeta < 0.012
					&& photonToTree.chargedIso < 2.6
					&& photonToTree.neutralIso < 3.5+0.04*photonToTree.pt
					&& photonToTree.photonIso < 1.3+0.005*photonToTree.pt
					&& photonToTree.neutralIso < 3.5+0.04*photonToTree.ptJet()
					&& photonToTree.photonIso < 1.3+0.005*photonToTree.ptJet()
				)
				// and the endcap definition, which is not used now
				|| ( std::abs( photonToTree.eta ) >= susy::etaGapEnd
					&& std::abs( photonToTree.eta ) <= susy::etaMax
					&& photonToTree.hadTowOverEm < 0.05
					&& photonToTree.sigmaIetaIeta < 0.034
					&& photonToTree.chargedIso < 2.3
					&& photonToTree.neutralIso < 2.9+0.04*photonToTree.pt
				);

			bool isPhoton = isPhotonOrElectron && !photonToTree.pixelseed;
			bool isPhotonElectron = isPhotonOrElectron && photonToTree.pixelseed;

			// photonJet definition
			bool isPhotonJet = !isPhotonOrElectron
				&& !photonToTree.pixelseed
				&& photonToTree.hadTowOverEm < 0.05
				&& photonToTree.sigmaIetaIeta < 0.012
				&& photonToTree.chargedIso < 26 && photonToTree.chargedIso > 0.26
				&& photonToTree.neutralIso < 35+0.4*photonToTree.pt && photonToTree.neutralIso > 0.35+0.004*photonToTree.pt
				&& photonToTree.photonIso < 13+0.05*photonToTree.pt && photonToTree.photonIso > 0.13+0.0005*photonToTree.pt
				&& (photonToTree.chargedIso < 5.2 || (photonToTree.neutralIso < 3.5 + 0.04*photonToTree.ptJet() && photonToTree.photonIso < 1.3 + 0.005*photonToTree.ptJet()))
				&& (photonToTree.neutralIso < 7 + 0.06*photonToTree.ptJet() || (photonToTree.chargedIso < 2.6 && photonToTree.photonIso < 1.3 + 0.005*photonToTree.ptJet()))
				&& (photonToTree.photonIso < 2.6 + 0.0075*photonToTree.ptJet() || (photonToTree.chargedIso < 2.6 && photonToTree.neutralIso < 3.5 + 0.04*photonToTree.ptJet()));

			// print photon information
			if( loggingVerbosity > 2 ) {
				if( isPhoton )         std::cout << " photon pT = " << photonToTree.pt << std::endl;
				if( isPhotonElectron ) std::cout << " photonElectron pT = " << photonToTree.pt << std::endl;
				if( isPhotonJet )      std::cout << " photonJet pT = " << photonToTree.pt << std::endl;
			}

			// Fill matching histograms only for photon-like objects
			if( !isPhotonOrElectron && !isPhotonJet ) continue;

			if( indexOfnearestParticle<tree::Particle>( photonToTree, genPhotons, .1, 0.9, 1.1, &hist2D["matchGenPhoton"] ) > -1 )
				photonToTree.setStatus( tree::kGenPhoton );
			if( indexOfnearestParticle<tree::Particle>( photonToTree, genElectrons, .1, -1e6, 1e6, &hist2D["matchGenPhoton"] ) > -1 )
				photonToTree.setStatus( tree::kGenElectron );

			// for plotting only
			const char* histname = "";
			if( isPhoton ) histname = "matchPhotonToJet";
			if( isPhotonJet ) histname = "matchPhotonJetToJet";
			if( isPhotonElectron ) histname = "matchPhotonElectronToJet";
			indexOfnearestParticle<tree::Jet>( photonToTree, jets, .2, .8, 3, &hist2D[histname] );

			if( loggingVerbosity > 2 )
				std::cout << "  ->jet pT = " << photonToTree._ptJet << std::endl;

			// If no jet is found for a loose photon, the photon is rejected
			// if( isPhotonJet && !photonToTree._ptJet ) continue;
			if( photonToTree.ptJet() < 110 ) continue;
			if( isPhoton )
				photons.push_back( photonToTree );
			if( isPhotonElectron )
				photonElectrons.push_back( photonToTree );
			if( isPhotonJet )
				photonJets.push_back( photonToTree );
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

		// this has to be done after the photon block, since leptons are cleared from photons
		fillLeptons();

		// do not allow leptons
		if( runType == kSimplifiedModel && ( electrons.size() || muons.size() ) ) continue;

		ht = getHt();
		nGoodJets = countGoodJets();
		if( loggingVerbosity > 1 ) {
			std::cout << "H_T = " << ht << std::endl;
			std::cout << "Found " << nGoodJets << " jets" << std::endl;
		}

		if( nGoodJets < 2 || ht < 500 ) continue;
		TVector3 mhtVector = getMhtVector();
		mht = mhtVector.Pt();
		mhtPhi = mhtVector.Phi();

		fillMetFilterBitHistogram( hist1D.at("metFilters"), event.metFilterBit );
		if( !event.passMetFilters() || !event.passMetFilter( susy::kEcalLaserCorr) ) continue;

		eventType eType = TreeWriter::whichEventType( photons, photonElectrons, photonJets );
		TVector3 recoilVector = getRecoilVector( eType );
		recoil = recoilVector.Pt();
		recoilPhi = recoilVector.Phi();

		if( eType == kPhotonEvent ) {
			photonTree.Fill();
			tryFill( hist1D, "gMet",signalPointString, met, weight );
			tryFill( hist1D, "gMetJesUp",signalPointString, met, weight );
			tryFill( hist1D, "gMetJesDown",signalPointString, met, weight );
			tryFill( hist1D, "gMetPuUp",signalPointString, met, weightPuUp );
			tryFill( hist1D, "gMetPuDown",signalPointString, met, weightPuDown );
			tryFill( hist1D, "gHt",signalPointString, ht, weight );
			tryFill( hist1D, "gNJets",signalPointString, nGoodJets, weight );
			tryFill( hist1D, "gPt",signalPointString, photons.at(0).ptJet(), weight );
		}
		if( eType == kJetEvent ) {
			photonJetTree.Fill();
			float qcdWeight=0, qcdWeightError=0;
			getQcdWeights( photonJets.at(0).ptJet(), recoil, qcdWeight, qcdWeightError );
			tryFill( hist1D, "fMet",signalPointString, met, weight*qcdWeight );
			tryFill( hist1D, "fMetError",signalPointString, met, weight*qcdWeightError );
		}
		if( eType == kElectronEvent ) {
			photonElectronTree.Fill();
			float ewkFakeRate = event.isRealData ?
				1. - 0.993 * (1. - std::pow(photonElectrons.at(0).pt / 2.9 + 1., -2.4)) * (1. - 0.23 * std::exp(-0.2777 * nTracksPV))* (1. - 5.66e-4 * nVertex)
				: 1 - (1 - 0.00623) * (1 - std::pow(photonElectrons.at(0).pt / 4.2 + 1,-2.9)) * (1 - 0.29 * std::exp(-0.335 * nTracksPV)) * (1 - 0.000223 * nVertex);
			tryFill( hist1D, "eMet",signalPointString, met, weight*ewkFakeRate );
		}

	} // for jentry

	outFile.cd();
	if( !jetScale ) {
		if ( runType == kTree || runType == kFullTree ) {
			photonTree.Write();
			photonElectronTree.Write();
			photonJetTree.Write();
			nPhotons.Write();
			eventNumbers.Write();
			for( std::map<std::string, TH2F>::const_iterator it = hist2D.begin();
					it!= hist2D.end(); ++it )
				it->second.Write();
		} // only for trees end

		// Write all histograms except for the Jec ones
		for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
					it!= hist1D.end(); ++it ) {
			if( ((std::string)(it->second.GetName())).find("Jes")==std::string::npos )
				it->second.Write();
		}
	} else if (jetScale == -1) {
		for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
					it!= hist1D.end(); ++it ) {
			if( ((std::string)(it->second.GetName())).find("JesDown")!=std::string::npos )
				it->second.Write();
		}
	}
	else if (jetScale == 1 ) {
		for( std::map<std::string, TH1F>::iterator it = hist1D.begin();
					it!= hist1D.end(); ++it ) {
			if( ((std::string)(it->second.GetName())).find("JesUp")!=std::string::npos )
				it->second.Write();
		}
	}

}


