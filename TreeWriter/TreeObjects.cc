#include "TreeObjects.h"

bool tree::Photon::isGen( genParticles id ) const {
	return genInformation & 1 << id;
}

void tree::Photon::setGen( genParticles id ) {
	genInformation |= 1 << id;
}

float tree::Photon::ptJet() const {
	if( _ptJet )
		return _ptJet;
	else
		return pt;
}

bool tree::EtGreater(const tree::Particle p1, const tree::Particle p2) {
	return p1.pt > p2.pt;
}

