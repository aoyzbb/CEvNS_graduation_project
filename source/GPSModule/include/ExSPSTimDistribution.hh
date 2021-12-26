//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        ExSPSTimDistribution.hh
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// 2016, Linyan, Created.
//    To enable time distribution storing and sampling
//
///////////////////////////////////////////////////////////////////////////////
//
//
// Class Description:
//
// To generate the time of a primary vertex according to the defined distribution 
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// ExSPSTimDistribution ()
//    Constructor: Initializes variables
//
// ~ExSPSTimDistribution ()
//    Destructor: 
//
// void SetTimeDisType(G4String)
//    Allows the user to choose the time distribution type. The arguments
//    are Mono (mono-energetic), Lin (linear), Pow (power-law), Exp 
//    (exponential), Gauss (gaussian), 
//    User (user-defined), Arb (arbitrary
//    point-wise), Epn (time per nucleon).
//
// void SetEmin(G4double)
//    Sets the minimum time.
//
// void SetEmax(G4double)
//    Sets the maximum time.
//
// void SetMonoTime(G4double)
//    Sets time for mono-energetic distribution.
//
// void SetAlpha(G4double)
//    Sets alpha for a power-law distribution.
//
//
// void SetEzero(G4double)
//    Sets Ezero for an exponential distribution.
//
// void SetGradient(G4double)
//    Sets gradient for a linear distribution.
//
// void SetInterCept(G4double)
//    Sets intercept for a linear distribution.
//
// void UserTimeHisto(G4ThreeVector)
//    Allows user to defined a histogram for the time distribution.
//
// void ArbTimeHisto(G4ThreeVector)
//    Allows the user to define an Arbitrary set of points for the
//    time distribution.
//
// void EpnTimeHisto(G4ThreeVector)
//    Allows the user to define an Time per nucleon histogram.
//
//
// void InputTimeSpectra(G4bool)
//    Allows the user to choose between momentum and time histograms
//    for user-defined histograms and arbitrary point-wise spectr.
//    The default is true (time).
//
// void InputDifferentialSpectra(G4bool)
//    Allows the user to choose between integral and differential 
//    distributions when using the arbitrary point-wise option.
//
// void ArbInterpolate(G4String)
//    ArbInterpolate allows the user to specify the type of function to
//    interpolate the Arbitrary points spectrum with.
//
//  void SetBiasRndm (ExSPSRandomGenerator* a)
//    Sets the biased random number generator
//
//  G4double GenerateOne(G4ParticleDefinition*);
//    Generate one random time for the specified particle
//
//  void ReSetHist(G4String);
//    Re-sets the histogram for user defined distribution
//
// void SetVerbosity(G4int)
//    Sets the verbosity level.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ExSPSTimDistribution_h
#define ExSPSTimDistribution_h 1

#include "G4PhysicsOrderedFreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "G4DataInterpolation.hh"

//
#include "G4SPSRandomGenerator.hh"

class ExSPSTimDistribution {
public:
	ExSPSTimDistribution();
	~ExSPSTimDistribution();

	static G4double fSegmentSize;

	G4bool bDependency;
	void SetTimeDisType(G4String);
	inline G4String GetTimeDisType() {
		return TimeDisType;
	}
	;
	void SetEmin(G4double);
	inline G4double GetEmin() {
		return Emin;
	}
	;
	inline G4double GetArbEmin() {
		return ArbEmin;
	}
	;
	void SetEmax(G4double);
	inline G4double GetEmax() {
		return Emax;
	}
	;
	inline G4double GetArbEmax() {
		return ArbEmax;
	}
	;
	void SetMonoTime(G4double);
	void SetAlpha(G4double);
	void SetBiasAlpha(G4double);
	void SetTemp(G4double);
	void SetBeamSigmaInE(G4double);
	void SetEzero(G4double);
	void SetGradient(G4double);
	void SetInterCept(G4double);
	void UserTimeHisto(G4ThreeVector);
	void ArbTimeHisto(G4ThreeVector);
	void ArbTimeHistoFile(G4String);
	void EpnTimeHisto(G4ThreeVector);

	void InputTimeSpectra(G4bool);
	void InputDifferentialSpectra(G4bool);
	void ArbInterpolate(G4String);
	inline G4String GetIntType() {
		return IntType;
	}
	;
	void Calculate();
	//
	void SetBiasRndm(G4SPSRandomGenerator* a) {
		eneRndm = a;
	}
	;
	// method to re-set the histograms
	void ReSetHist(G4String);
	// Set the verbosity level.
	void SetVerbosity(G4int a) {
		verbosityLevel = a;
	}
	;
	//x
	G4double GetWeight() {
		return weight;
	}

	G4double GetMonoTime() {
		return MonoTime;
	}
	; //Mono-energteic time
	G4double GetSE() {
		return SE;
	}
	; // Standard deviation for Gaussion distrbution in time
	G4double Getalpha() {
		return alpha;
	}
	; // alpha (pow)
	G4double GetEzero() {
		return Ezero;
	}
	; // E0 (exp)
	G4double Getgrad() {
		return grad;
	}
	; // gradient and intercept for linear spectra
	G4double Getcept() {
		return cept;
	}
	;

	inline G4PhysicsOrderedFreeVector GetUserDefinedTimeHisto() {
		return UDefTimeH;
	}
	;
	inline G4PhysicsOrderedFreeVector GetArbTimeHisto() {
		return ArbTimeH;
	}
	;

	G4double GenerateOne(G4ParticleDefinition*);
	G4double GetProbability (G4double);


private:
	void LinearInterpolation();
	void LogInterpolation();
	void ExpInterpolation();
	void SplineInterpolation();

	// The following methods generate energies according to the spectral
	// parameters defined above.
	void GenerateMonoTime();
	void GenerateLinearTimes(G4bool);
	void GeneratePowTimes(G4bool);
	void GenerateBiasPowTimes();
	void GenerateExpTimes(G4bool);
	void GenerateGaussTimes();
	void GenUserHistTimes();
	void GenEpnHistTimes();
	void GenArbPointTimes();
	// converts time per nucleon to time.
	void ConvertEPNToTime();


private:

	G4String TimeDisType; // time dis type Variable  - Mono,Lin,Exp,etc
	G4double weight; // particle weight
	G4double MonoTime; //Mono-energteic time
	G4double SE; // Standard deviation for Gaussion distrbution in time
	G4double Emin, Emax; // emin and emax
	G4double alpha, Ezero; // alpha (pow), E0 (exp) 
	G4double biasalpha; // biased power index
	G4double grad, cept; // gradient and intercept for linear spectra
        G4double prob_norm; // normalisation factor use in calculate the probability 
        G4bool Biased; // true - biased to power-law
	G4bool TimeSpec; // true - time spectra, false - momentum spectra
	G4bool DiffSpec; // true - differential spec, false integral spec
	//G4bool ApplyRig; // false no rigidity cutoff, true then apply one
	//G4double ERig; // time of rigidity cutoff
	G4PhysicsOrderedFreeVector UDefTimeH; // time hist data
	G4PhysicsOrderedFreeVector IPDFTimeH;
	G4bool IPDFTimeExist, IPDFArbExist, Epnflag;
	G4PhysicsOrderedFreeVector ArbTimeH; // Arb x,y histogram
	G4PhysicsOrderedFreeVector IPDFArbTimeH; // IPDF for Arb
	G4PhysicsOrderedFreeVector EpnTimeH;
	G4double BBHist[10001];
	G4String IntType; // Interpolation type
	G4double Arb_grad[1024], Arb_cept[1024]; // grad and cept for 1024 segments
	G4double Arb_alpha[1024], Arb_Const[1024]; // alpha and constants
	G4double Arb_ezero[1024]; // ezero
	G4double ArbEmin, ArbEmax; // Emin and Emax for the whole arb distribution used primarily for debug.

	G4double particle_time;
	G4ParticleDefinition* particle_definition;

	G4SPSRandomGenerator* eneRndm;

	// Verbosity
	G4int verbosityLevel;

	G4PhysicsOrderedFreeVector ZeroPhysVector; // for re-set only

	G4DataInterpolation *SplineInt[1024]; // holds Spline stuff required for sampling
	G4DataInterpolation *Splinetemp; // holds a temp Spline used for calculating area

};

#endif

