//
// Improved from BetaSpec.
// By Guo Ziyi 2017.5
//
#ifndef ALPHABRANCH_H
#define ALPHABRANCH_H

#include <string>
#include <vector>
#include "DecayBranch.hh"


class AlphaBranch : public DecayBranch
{
public:
    AlphaBranch();
	AlphaBranch(double ratio, double newEk, std::vector<double>& gammalist);

    ~AlphaBranch();


	// Get random Ek from spectrum, in MeV
	double GetRandomEk();

private:
	double Ek;

};


#endif
