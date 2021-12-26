#include "AlphaBranch.hh"
using namespace std;


AlphaBranch::AlphaBranch()
{
	Ek = 0;
	type = Alpha;
	fBranchRatio = 0;

}

AlphaBranch::AlphaBranch(double ratio, double newEk, vector<double>& gammalist)
{
	type = Alpha;
	Ek = newEk;
	fBranchRatio = ratio;
	fGammaEnergy = gammalist;

}

AlphaBranch::~AlphaBranch()
{

}

double AlphaBranch::GetRandomEk()
{
	return Ek;
}
