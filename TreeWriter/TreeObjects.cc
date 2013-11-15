#include "TreeObjects.h"

bool tree::Photon::isGen( genParticles id ) const {
	return genInformation & 1 << id;
}

void tree::Photon::setGen( genParticles id ) {
	genInformation |= 1 << id;
}

bool tree::Jet::isMatch( jetMatches id ) const {
	return matchInformation & 1 << id;
}

void tree::Jet::setMatch( jetMatches id ) {
	matchInformation |=1 << id;
}

float tree::Photon::ptJet() const {
	return _ptJet ? _ptJet : pt;
}

bool tree::EtGreater(const tree::Particle p1, const tree::Particle p2) {
	return p1.pt > p2.pt;
}

float tree::Particle::DeltaR( const Particle &p2 ) const {
	TVector2 tempVect;
	float dphi = tempVect.Phi_mpi_pi( phi -p2.phi );
	float deta = eta - p2.eta;
	return sqrt( dphi*dphi + deta*deta );
}

float tree::Particle::DeltaR( const TLorentzVector &vec2 ) const {
	TVector3 vec1;
	vec1.SetPtEtaPhi( 1, eta, phi );
	return vec1.DeltaR( vec2.Vect() );
}
