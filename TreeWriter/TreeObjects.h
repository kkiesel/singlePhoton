#ifndef TreeObjects_h
#define TreeObjects_h

namespace tree {
// In this namespace classes for the trees are defined.

class Particle {
	public:
		float pt, eta, phi;
};

class Photon : public Particle {
	public:
		float r9, sigmaIetaIeta, hadTowOverEm;
		float chargedIso, neutralIso, photonIso;
		bool conversionSafeVeto;
		int pixelseed;
};

class Jet : public Particle{
	public:
		float bCSV;
		float chargedHadronEnergy,neutralHadronEnergy,photonEnergy,electronEnergy,muonEnergy,HFHadronEnergy,HFEMEnergy,chargedEmEnergy,chargedMuEnergy,neutralEmEnergy;
};

bool EtGreater(const tree::Particle, const tree::Particle);

} // end namespace definition

#endif



