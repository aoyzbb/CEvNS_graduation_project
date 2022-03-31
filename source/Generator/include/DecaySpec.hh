//
// Improved from BetaSpec.
// By Guo Ziyi 2017.5
//
#ifndef DECAYSPEC_H
#define DECAYSPEC_H

#include <string>
#include <vector>
#include "DecayBranch.hh"

class TH1D;

class DecaySpec {
public:
    DecaySpec();
    DecaySpec(std::string name, int z, int a);
    virtual ~DecaySpec();


	DecayType GetVertex(double& Ek, std::vector<double>& eGammas);

	// Use for root idle
	TH1D* GetAlphaEkSpe();
	TH1D* GetAlphaEkSpeWithGamma();
	TH1D* GetBetaMinusEkSpe();
	TH1D* GetBetaMinusEkSpeWithGamma();
    
	std::vector<DecayBranch*> DecayBranchList;

private:
	std::string isotope;
	int Z;
	int A;
	double m_BrNom;

protected:
	void LoadDataFile(std::string filename, DecayType type);

};

#endif
