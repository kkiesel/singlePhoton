#ifndef TreeObjects_h
#define TreeObjects_h

namespace tree {
// In this namespace classes for the trees are defined.

enum genParticles{
	genPhoton,
	genElectron
};

class Particle {
	public:
		float pt, eta, phi;
};

class Photon : public Particle {
	public:
		bool isGen( genParticles id ) const;
		void setGen( genParticles id );
		float ptJet() const;
		float _ptJet;
		float sigmaIphiIphi;
		float r9, sigmaIetaIeta, hadTowOverEm;
		float chargedIso, neutralIso, photonIso;
		bool conversionSafeVeto;
		int pixelseed;
		short genInformation;
};

class Jet : public Particle{
	public:
		float bCSV;
		float chargedHadronEnergy,neutralHadronEnergy,photonEnergy,electronEnergy,muonEnergy,HFHadronEnergy,HFEMEnergy,chargedEmEnergy,chargedMuEnergy,neutralEmEnergy;
};

bool EtGreater(const tree::Particle, const tree::Particle);

} // end namespace definition

#endif



