#include "treeWriter.h"

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

bool isLooseJet( const susy::PFJet& jet ) {
  /**
   * \brief Apply loose cut on jets.
   *
   * See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_7_TeV_data_a
   * for more information.
   */
  double energy = jet.momentum.E();
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

///////////////////////////////////////////////////////////////////////////////
// Here the class implementation begins ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

TreeWriter::TreeWriter( int nFiles, char** fileList, std::string const& outputName ) :
  inputTree("susyTree"),
  event(),
  outFile( outputName.c_str(), "recreate" ),
  outputTree("myTree","Tree for single photon analysis")
{

  for( int i = 0; i<nFiles; ++i )
    inputTree.Add( fileList[i] );

  event.setInput( inputTree );

  outputTree.Branch("jets", &jets);
  outputTree.Branch("met", &met, "met/F");

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
      break;
    }
  }

  if( loggingVerbosity > 2 )
    std::cout << "Pile-up weight = " << thisWeight << std::endl;
  return thisWeight;
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
    if( corrP4.Pt() < 40 ) continue;
    if( !isLooseJet( *it ) ) continue;

    jetToTree.pt = corrP4.Pt();
    jetToTree.eta = corrP4.Eta();
    jetToTree.phi = corrP4.Phi();
    jets.push_back( jetToTree );

    if( loggingVerbosity > 2 )
      std::cout << " p_T, jet = " << jetToTree.pt << std::endl;
  }// for jet

  if( loggingVerbosity > 1 )
      std::cout << "Found " << jets.size() << " uncleaned jets" << std::endl;
}

void TreeWriter::Loop() {
  /**
   * \brief Loops over input chain and fills tree
   *
   * This is the major function of treeWriter, which initialize the output, loops
   * over all events and fill the tree. In the end, the tree is saved to the
   * output File
   */

  for (long jentry=0; jentry < inputTree.GetEntries(); ++jentry) {
    if( jentry > 100 ) break;
    event.getEntry(jentry);

    // For data, the weight is 1. Else take the pileup weight.
    weight = event.isRealData ? 1. : getPileUpWeight();

    if ( event.isRealData && !isGoodLumi() ) continue;
    if ( event.isRealData && !passTrigger() ) continue;

    susy::MET pfMet = event.metMap["pfMet"];
    met = pfMet.met();

    fillJets();

    outputTree.Fill();
  }


  outFile.cd();
  outputTree.Write();
}


