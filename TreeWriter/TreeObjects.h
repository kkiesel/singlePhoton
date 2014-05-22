#include <TLorentzVector.h>
#include <TVector3.h>

#ifndef TreeObjects_h
#define TreeObjects_h

namespace tree {
// In this namespace classes for the trees are defined.

enum genParticles {
	kGenPhoton,
	kGenElectron,
	kGenJet,
	kNearLepton
};

enum jetMatches {
	kJetId,
	kJetPhoton,
	kJetAllPhoton,
	kJetCount
};

class Particle {
	public:
		float DeltaR( const Particle &p2 ) const;
		float DeltaR( const TLorentzVector &vec2 ) const;
		float DeltaPhi( float phi2 ) const;
		float pt, eta, phi;
};

class Photon : public Particle {
	public:
		bool isGen( genParticles id ) const;
		void setGen( genParticles id );
		float ptJet() const;
		float _ptJet;
		float _etaJet;
		float _phiJet;
		float sigmaIphiIphi;
		float r9, sigmaIetaIeta, hadTowOverEm;
		float chargedIso, neutralIso, photonIso;
		bool conversionSafeVeto;
		int pixelseed;
		short genInformation;
		short matchedJetIndex;
};

class Jet : public Particle{
	public:
		//float DeltaR( const Photon &p2 ) const;
		//float DeltaR( const Jet &p2 ) const;
		bool isMatch( jetMatches id ) const;
		void setMatch( jetMatches id );
		short matchInformation;
		float bCSV;
		float chargedHadronEnergy,
			neutralHadronEnergy,
			photonEnergy,
			electronEnergy,
			muonEnergy,
			HFHadronEnergy,
			HFEMEnergy,
			chargedEmEnergy,
			chargedMuEnergy,
			neutralEmEnergy;
};

bool EtGreater(const tree::Particle, const tree::Particle);

} // end namespace definition

#endif



