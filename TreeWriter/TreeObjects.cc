#include "TreeObjects.h"

float tree::Photon::ptJet() const {
	if( _ptJet )
		return _ptJet;
	else
		return pt;
}

bool tree::EtGreater(const tree::Particle p1, const tree::Particle p2) {
	return p1.pt > p2.pt;
}

