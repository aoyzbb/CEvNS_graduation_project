//
// Improved from BetaSpec.
// By Guo Ziyi 2017.5
//
#ifndef DECAYBRANCH_H
#define DECAYBRANCH_H

#include <string>
#include <vector>
#include "DecayBranch.hh"


enum DecayType
{
	Alpha = 0,
	BetaMinus,
	BetaPlus,
	EC,
	Gamma,
	INVALID
};

class DecayBranch 
{
public:
    DecayBranch();
    virtual ~DecayBranch();

	// Get decay type
	DecayType GetType() { return type; }
	// Get branch ratio
	double GetBranchRatio() { return fBranchRatio; }

	// Get gamma energy, in MeV
	std::vector<double> fGammaEnergy;
	// Get random Ek from spectrum, in MeV
	virtual double GetRandomEk();

	static const char* DecayTypeString[];

protected:
	DecayType type;
	double fBranchRatio;



};

#endif
