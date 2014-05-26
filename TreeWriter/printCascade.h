void printChildren( int index, const susy::ParticleCollection&  particles, int level=0 ){
	/* Prints the daughter and recursivly her's daughters of particle 'index'.
	 * The status of the particle is denoted in paranthesis.
	 */
	for (int i = 0; i< level; ++i )
		std::cout <<"\t";

	std::map<int,std::string> pdgIdToString;
	std::map<int,std::string> pdgIdToStringTmp;
	pdgIdToString[1] = "d";
	pdgIdToString[2] = "u";
	pdgIdToString[3] = "s";
	pdgIdToString[4] = "c";
	pdgIdToString[5] = "b";
	pdgIdToString[6] = "t";
	pdgIdToString[11] = "e-";
	pdgIdToString[12] = "νe_-";
	pdgIdToString[13] = "μ-";
	pdgIdToString[14] = "νμ-";
	pdgIdToString[15] = "τ-";
	pdgIdToString[16] = "ντ-";
	pdgIdToString[21] = "g";
	pdgIdToString[22] = "γ";
	pdgIdToString[23] = "Z";
	pdgIdToString[24] = "W";
	pdgIdToString[221] = "η";
	pdgIdToString[111] = "π";
	pdgIdToString[211] = "π+";
	pdgIdToString[2212] = "p";

	for( std::map<int,std::string>::iterator it = pdgIdToString.begin(); it != pdgIdToString.end(); ++it ) {
		pdgIdToStringTmp[it->first] = it->second;
		pdgIdToStringTmp[1000000+it->first] = std::string("~")+it->second;
	}
	pdgIdToStringTmp[1000022] = "~χ1";
	pdgIdToStringTmp[1000023] = "~χ2";
	pdgIdToStringTmp[1000025] = "~χ3";
	pdgIdToStringTmp[1000035] = "~χ4";
	pdgIdToStringTmp[1000039] = "~G";
	pdgIdToString = pdgIdToStringTmp;
	for( std::map<int,std::string>::iterator it = pdgIdToString.begin(); it != pdgIdToString.end(); ++it ) {
		pdgIdToStringTmp[it->first] = it->second;
		pdgIdToStringTmp[-it->first] = it->second + std::string("-");
	}
	pdgIdToStringTmp[1000024] = "~χ1+";
	pdgIdToStringTmp[-1000024] = "~χ1-";
	pdgIdToStringTmp[1000037] = "~χ2+";
	pdgIdToStringTmp[-1000037] = "~χ2-";
	pdgIdToString = pdgIdToStringTmp;

	if( pdgIdToString.find( particles[index].pdgId) != pdgIdToString.end() )
		std::cout << pdgIdToString.find( particles[index].pdgId)->second;
	else
		std::cout << particles[index].pdgId;

	std::cout << " (" << (int)particles[index].status << ") " << particles[index].momentum.Pt() << std::endl;

	//susy::ParticleCollection children;
	for( susy::ParticleCollection::const_iterator it = particles.begin(); it != particles.end(); ++it ) {
		if( it->motherIndex == index )
			printChildren( std::distance<susy::ParticleCollection::const_iterator>(particles.begin(), it ), particles, level+1 );
	}
}

void printCascade( const susy::ParticleCollection & genParticles ) {

	std::cout << "=================================================" << std::endl;

	for( susy::ParticleCollection::const_iterator it = genParticles.begin();
			it != genParticles.end(); ++it ){
		if( it->motherIndex == -1 ){
			printChildren( std::distance<susy::ParticleCollection::const_iterator>(genParticles.begin(), it ), genParticles );
		}
	}
}

