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
		float r9, sigmaIetaIeta, hadTowOverEm, pixelseed;
		float chargedIso, neutralIso, photonIso;
};

class Jet : public Particle{
	public:
		float pt, eta;
};

bool EtGreater(const tree::Particle, const tree::Particle);


} // end namespace definition

#endif
