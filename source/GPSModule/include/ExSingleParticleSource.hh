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
// MODULE:        ExSingleParticleSource.hh
//
// Version:      1.0
// Date:         5/02/04
// Author:       Fan Lei
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//      Based on the ExGeneralParticleSource class in Geant4 v6.0
//
// 2016, Linyan, Modified
// To support predefine dradioactivity and timestamp
//
///////////////////////////////////////////////////////////////////////////////
//
// Class Description:
//
// The Single Particle Source is designed to extend the functionality of the
// G4ParticleGun class. It is designed to allow specification of input
// particles in terms of position, direction (or angular) and energy
// distributions.  It is used by the General Particle source class
// and it is derived from G4VPrimaryGenerator.
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// ExSingleParticleSource ()
//    Constructor: Initializes variables and instantiates the
//                 Messenger and Navigator classes
//
// ~ExSingleParticleSource ()
//    Destructor:  deletes Messenger and prints out run information.
//
// void GeneratePrimaryVertex(G4Event *evt)
//    Generate the particles initial parameters.
//
//  ExSPSPosDistribution* GetPosDist()
//    Return a pointer to the position distribution generator
//
//  ExSPSAngDistribution* GetAngDist()
//    Return a pointer to the angular distribution generator
//
//  G4SPSEneDistribution* GetEneDist()
//     Return a pointer to the energy distribution generator
//
//  G4SPSRandomGenerator* GetBiasRndm() {return biasRndm;};
//     Return a pointer to the biased random number generator
//
//  void SetVerbosity(G4int);
//     Set the verbosity level.
//
//  void SetParticleDefinition ();
//  G4ParticleDefinition * GetParticleDefinition ()
//     Get/Set the particle definition of the primary track
//
//  void SetParticleCharge(G4double aCharge)
//     set the charge state of the primary track
//
//  inline void SetParticlePolarization (G4ThreeVector aVal)
//  inline G4ThreeVector GetParticlePolarization ()
//     Set/Get the polarization state of the primary track
//
//  inline void SetParticleTime(G4double aTime)  { particle_time = aTime; };
//  inline G4double GetParticleTime()  { return particle_time; };
//     Set/Get the Time.
//
//  inline void SetNumberOfParticles(G4int i)
//  inline G4int GetNumberOfParticles()
//     set/get the number of particles to be generated in the primary track
//
//  inline G4ThreeVector GetParticlePosition()
//  inline G4ThreeVector GetParticleMomentumDirection()
//  inline G4double GetParticleEnergy()
//     get the position, direction, and energy of the current particle
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef ExSingleParticleSource_h
#define ExSingleParticleSource_h 1

#include "globals.hh"
#include <vector>

//#include "TimeStamp.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
//
#include "ExSPSPosDistribution.hh"
#include "ExSPSAngDistribution.hh"
#include "ExSPSEneDistribution.hh"
#include "ExSPSTimDistribution.hh"
#include "G4SPSRandomGenerator.hh"
#include "ExGeneralParticleSourceMessenger.hh"
#include "G4Geantino.hh"

#include "BetaSpec.hh"
#include "DecaySpec.hh"

struct SingleParticle
{
	G4String fName;
	int RadZ;
	int RadA;
	G4ParticleDefinition *fDef;
	DecaySpec *fDecaySpec;
	G4bool UserRad;
	G4ParticleMomentum fMomentumDirection;
	G4double fEnergy;
	G4double fCharge;
	G4ThreeVector fPosition;
	G4ThreeVector fPolarization;
	G4double fWeight;
	//TimeStamp fTimeStamp;

	ExSPSPosDistribution *posGenerator;
	ExSPSAngDistribution *angGenerator;
	ExSPSEneDistribution *eneGenerator;
	ExSPSTimDistribution *timGenerator;
	ExSPSTimDistribution *timsecGenerator;
	ExSPSTimDistribution *timnsecGenerator;

	SingleParticle(G4SPSRandomGenerator *biasRndm)
	{
		fName = "geantino";
		UserRad = false;
		RadZ = 0;
		RadA = 0;
		fDef = G4Geantino::Geantino();
		fMomentumDirection = G4ParticleMomentum(1, 0, 0);
		fEnergy = 0;
		fCharge = 0.0;
		G4ThreeVector zero;
		fPosition = zero;
		fPolarization = zero;
		//fTimeStamp = TimeStamp(0);
		fDecaySpec = new DecaySpec;
		posGenerator = new ExSPSPosDistribution();
		posGenerator->SetBiasRndm(biasRndm);
		angGenerator = new ExSPSAngDistribution();
		angGenerator->SetBiasRndm(biasRndm);
		angGenerator->SetPosDistribution(posGenerator);
		angGenerator->SetAngDistType("iso");

		eneGenerator = new ExSPSEneDistribution();
		eneGenerator->SetBiasRndm(biasRndm);
		timGenerator = new ExSPSTimDistribution();
		timGenerator->SetBiasRndm(biasRndm);
		timsecGenerator = new ExSPSTimDistribution();
		timsecGenerator->SetBiasRndm(biasRndm);
		timnsecGenerator = new ExSPSTimDistribution();
		timnsecGenerator->SetBiasRndm(biasRndm);
	}

	~SingleParticle()
	{
	}
};

class ExSingleParticleSource : public G4VPrimaryGenerator
{
public:
	ExSingleParticleSource(ExGeneralParticleSourceMessenger *theM);
	~ExSingleParticleSource();
	void GeneratePrimaryVertex(G4Event *evt);
	//

	G4int GetCurrentSourceIndex() { return currentIdx; };
	G4int GetCurrentParticleIndex() { return currentParticleIdx; };
	ExSPSPosDistribution *GetPosDist()
	{
		return posGenerator[currentIdx];
	}
	ExSPSPosDistribution *GetCurrentPosDist()
	{
		return particleList[currentIdx]->posGenerator;
	}

	ExSPSAngDistribution *GetAngDist()
	{
		return angGenerator[currentIdx];
	}

	ExSPSAngDistribution *GetCurrentAngDist()
	{
		return particleList[currentIdx]->angGenerator;
	}

	ExSPSEneDistribution *GetCurrentEneDist()
	{
		return particleList[currentIdx]->eneGenerator;
	}

	ExSPSTimDistribution *GetCurrentTimDist()
	{
		return particleList[currentIdx]->timGenerator;
	}

	ExSPSTimDistribution *GetCurrentSecTimDist()
	{
		return particleList[currentIdx]->timsecGenerator;
	}

	ExSPSTimDistribution *GetCurrentNSecTimDist()
	{
		return particleList[currentIdx]->timnsecGenerator;
	}

	G4SPSRandomGenerator *GetBiasRndm()
	{
		return biasRndm;
	};

	//void ClearTimeStampList()
	//{
	//	fTotalTimeStamp.clear();
	//}

	G4int GetTotalIdx() { return NumberOfParticlesToBeGenerated.size(); };

	//G4int GetTotalParticleNumber() { return fTotalTimeStamp.size(); }
	// Set the verbosity level.
	void SetVerbosity(G4int);
	G4int GetVerbosityLevel() { return verbosityLevel; };

	// Set the particle species
	void SetParticleDefinition(G4ParticleDefinition *aParticleDefinition);
	void SetParticleDefinition(G4String name);
	void SetParticleDefinition(int Z, int A);
	// void SetNewPrimary(G4PrimaryParticle* particle, G4int index);
	void SetCurrentNewPrimary(G4PrimaryParticle *particle, G4int index, G4ParticleDefinition *fDef, G4ThreeVector fPosition, G4double fEnergy, G4ParticleMomentum fMomentumDirection, G4double fCharge, G4ThreeVector fPolarization);
	void SetCurrentNewPrimary(G4PrimaryParticle *particle, G4int index);
	void SetRadZA(int Z, int A);
	inline G4ParticleDefinition *GetParticleDefinition()
	{
		return particle_definition[currentIdx];
	};

	inline void SetParticleCharge(G4double aCharge)
	{
		particle_charge[currentIdx] = aCharge;
	};

	// Set polarization
	inline void SetParticlePolarization(G4ThreeVector aVal)
	{
		particle_polarization[currentIdx] = aVal;
	};
	inline G4ThreeVector GetParticlePolarization()
	{
		return particle_polarization[currentIdx];
	};

	// Set Time.
	//inline void SetParticleTimeStamp(G4int sec, G4int msec)
	//{
	//	TimeStamp newstamp(sec);
	//	newstamp.Add(msec * 1e-9);
	//	particle_stamp.push_back(newstamp);
	//}
	//inline TimeStamp GetParticleTimeStamp(G4int nI)
	//{
	//	return particle_stamp[nI];
	//}
	//inline TimeStamp GetCurrentParticleTimeStamp(G4int nI)
	//{
	//	return fTotalTimeStamp[nI];
	//}

	//;

	inline void SetNumberOfParticles(G4int i)
	{
		NumberOfParticlesToBeGenerated[currentIdx] = i;
		nParticlesToBeGen = i;
	};
	//
	inline G4int GetNumberOfParticles()
	{
		return NumberOfParticlesToBeGenerated[currentIdx];
	};
	inline G4ThreeVector GetParticlePosition()
	{
		return particle_position[currentIdx];
	};
	inline G4ThreeVector GetParticleMomentumDirection()
	{
		return particle_momentum_direction[currentIdx];
	};
	inline G4double GetParticleEnergy()
	{
		return particle_energy[currentIdx];
	};
	void AddaEvent(G4double aV);

	//	void CopyPosDistribution(ExSPSPosDistribution* pos0, ExSPSPosDistribution* pos1);

private:
	//std::vector<TimeStamp> particle_stamp;

	std::vector<ExSPSPosDistribution *> posGenerator;
	std::vector<ExSPSAngDistribution *> angGenerator;
	std::vector<ExSPSEneDistribution *> eneGenerator;
	std::vector<ExSPSTimDistribution *> timGenerator;
	std::vector<ExSPSTimDistribution *> timsecGenerator;
	std::vector<ExSPSTimDistribution *> timnsecGenerator;
	G4SPSRandomGenerator *biasRndm;
	//
	// Other particle properties
	std::vector<G4int> NumberOfParticlesToBeGenerated;
	std::vector<G4ParticleDefinition *> particle_definition;
	std::vector<G4ParticleMomentum> particle_momentum_direction;
	std::vector<G4double> particle_energy;
	std::vector<G4double> particle_charge;
	std::vector<G4ThreeVector> particle_position;
	std::vector<G4ThreeVector> particle_polarization;
	std::vector<G4double> particle_weight;

	G4int currentIdx;
	ExGeneralParticleSourceMessenger *theMessenger;
	// Verbosity
	G4int verbosityLevel;

	// New style
	G4int nParticlesToBeGen;
	G4int currentParticleIdx;
	std::vector<SingleParticle *> particleList;
	//std::vector<TimeStamp> fTotalTimeStamp;

	G4int RadZ, RadA;

	std::vector<SingleParticle> fParticleList;
};

#endif
