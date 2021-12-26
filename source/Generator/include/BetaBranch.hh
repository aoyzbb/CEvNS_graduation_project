//
// Improved from BetaSpec.
// By Guo Ziyi 2017.5
//
#ifndef BETABRANCH_H
#define BETABRANCH_H

#include <string>
#include <vector>
#include "DecayBranch.hh"
#include "TComplex.h"

class TGraph;
class TH1D;

class BetaBranch : public DecayBranch
{
public:
    BetaBranch();
	BetaBranch(int branchId, double ratio, DecayType newtype, double Qv, int z, int a, int ds, int dp, std::vector<double>& gammalist);

    ~BetaBranch();


	// Get random Ek from spectrum, in MeV
	double GetRandomEk();

private:
	double Qvalue;
	int Z;
	int A;
	int C;
	int Zdaughter;
	int dSpin;
	int dParity;
    //double Qmax;    // maximum energy level difference
    //double Emax;    // maximum (energy level difference + all gamma energies)
    double R;       // nucleus charge radius
    int verbosity;

    static double me;		// mass of electron
    static double alpha;	// fine structure constant
    TGraph* gScreeningV;
	TH1D* hBeta;
	double EnergyHisMax;
	//std::vector<double> fSpectraVector;
    bool useForbiddennessHuber;
	
protected:
	static TComplex CGamma(TComplex z);	// Complex Gamma function

    //static int WeightedChoice(double* weights, int size);

	// correction functions
	bool UseForbiddenCorrection(int deltaSpin, int deltaParity);
	double Fermi(int z, double Te, int charge=-1);
	double Screening(int z, double Te, int charge=-1);
	double Forbiddenness(int z, double Te, double Q0, int deltaSpin, int deltaParity, int charge=-1);
	double FiniteSizeEM(int z, double Te, int charge=-1);
    double FiniteSizeWI(int z, double Te, double Q0, int charge=-1);
    double WeakMagnetism(double Te, int charge=-1);



};


#endif
