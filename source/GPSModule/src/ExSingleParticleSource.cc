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
#include <cmath>
#include <vector>

#include "ExSingleParticleSource.hh"
#include "ExGeneralParticleSourceMessenger.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "JPSimPrimaryVertexInformation.hh"
#include "TH1.h"
//#include "TimeStamp.hh"

ExSingleParticleSource::ExSingleParticleSource(ExGeneralParticleSourceMessenger *theM)
{
	// Initialise all variables
	// Position distribution Variables
	theMessenger = theM;

	NumberOfParticlesToBeGenerated.clear();
	particle_definition.clear();
	particle_momentum_direction.clear();
	particle_energy.clear();
	particle_position.clear();
	particle_polarization.clear();
	particle_charge.clear();
	particle_weight.clear();
	posGenerator.clear();
	angGenerator.clear();
	eneGenerator.clear();
	timGenerator.clear();
	timsecGenerator.clear();
	timnsecGenerator.clear();
	currentIdx = 0;

	biasRndm = new G4SPSRandomGenerator();

	particleList.push_back(new SingleParticle(biasRndm));
	nParticlesToBeGen = 1;
	currentParticleIdx = 0;

	NumberOfParticlesToBeGenerated.push_back(1);
	particle_definition.push_back(G4Geantino::GeantinoDefinition());
	G4ThreeVector zero;
	particle_momentum_direction.push_back(G4ParticleMomentum(1, 0, 0));
	particle_energy.push_back(1.0 * MeV);
	particle_position.push_back(zero);
	particle_polarization.push_back(zero);
	particle_charge.push_back(0.0);
	particle_weight.push_back(1.0);
	//particle_stamp.push_back(TimeStamp(0));

	posGenerator.push_back(new ExSPSPosDistribution());
	posGenerator[currentIdx]->SetBiasRndm(biasRndm);
	angGenerator.push_back(new ExSPSAngDistribution());
	angGenerator[currentIdx]->SetPosDistribution(posGenerator[currentIdx]);
	angGenerator[currentIdx]->SetBiasRndm(biasRndm);
	eneGenerator.push_back(new ExSPSEneDistribution());
	eneGenerator[currentIdx]->SetBiasRndm(biasRndm);
	timGenerator.push_back(new ExSPSTimDistribution());
	timGenerator[currentIdx]->SetBiasRndm(biasRndm);
	timsecGenerator.push_back(new ExSPSTimDistribution());
	timsecGenerator[currentIdx]->SetBiasRndm(biasRndm);
	timnsecGenerator.push_back(new ExSPSTimDistribution());
	timnsecGenerator[currentIdx]->SetBiasRndm(biasRndm);

	currentIdx = G4int(NumberOfParticlesToBeGenerated.size() - 1);
	// verbosity
	verbosityLevel = 0;
}

ExSingleParticleSource::~ExSingleParticleSource()
{
	delete biasRndm;
}

void ExSingleParticleSource::AddaEvent(G4double aV)
{
	if (verbosityLevel >= 1)
		G4cout << "add correlated event: " << aV << G4endl;

	// G4cout<<aV<<G4endl;
	NumberOfParticlesToBeGenerated.push_back(int(aV));
	nParticlesToBeGen++;
	currentIdx = G4int(NumberOfParticlesToBeGenerated.size() - 1);
	currentParticleIdx = nParticlesToBeGen - 1;
	// G4cout<<"AddA: "<<currentIdx<<G4endl;

	// G4cout<<"AddaEvent"<<G4endl;
	particleList.push_back(new SingleParticle(biasRndm));
	particleList[currentParticleIdx]->posGenerator->CopyFromSPSP(particleList[currentParticleIdx - 1]->posGenerator);
	particleList[currentParticleIdx]->angGenerator->SetAngDistType("iso");

	// G4cout<<"AddaEventEnd"<<G4endl;

	particle_definition.push_back(G4Geantino::GeantinoDefinition());
	G4ThreeVector zero;
	particle_momentum_direction.push_back(G4ParticleMomentum(1, 0, 0));
	particle_energy.push_back(1.0 * MeV);
	particle_position.push_back(zero);
	particle_polarization.push_back(zero);
	particle_charge.push_back(0.0);
	particle_weight.push_back(1.0);
	//particle_stamp.push_back(TimeStamp(0));

	// G4cout<<"AddA posGenerator size "<<posGenerator.size()<<G4endl;
	posGenerator.push_back(new ExSPSPosDistribution());
	// G4cout<<"AddA posGenerator after size "<<posGenerator.size()<<G4endl;
	posGenerator[currentIdx]->SetBiasRndm(biasRndm);
	// G4cout<<"AddA currentIdx "<<currentIdx<<G4endl;
	// G4cout<<"AddA Last Confine: "<<currentIdx-1<<" "<<posGenerator[currentIdx-1]->GetConfineVolume()<<G4endl;
	posGenerator[currentIdx]->CopyFromSPSP(posGenerator[currentIdx - 1]);

	angGenerator.push_back(new ExSPSAngDistribution());
	angGenerator[currentIdx]->SetPosDistribution(posGenerator[currentIdx]);
	angGenerator[currentIdx]->SetAngDistType("iso");
	angGenerator[currentIdx]->SetBiasRndm(biasRndm);
	eneGenerator.push_back(new ExSPSEneDistribution());
	eneGenerator[currentIdx]->SetBiasRndm(biasRndm);
	timGenerator.push_back(new ExSPSTimDistribution());
	timGenerator[currentIdx]->SetBiasRndm(biasRndm);
	timsecGenerator.push_back(new ExSPSTimDistribution());
	timsecGenerator[currentIdx]->SetBiasRndm(biasRndm);
	timnsecGenerator.push_back(new ExSPSTimDistribution());
	timnsecGenerator[currentIdx]->SetBiasRndm(biasRndm);

	if (verbosityLevel >= 1)
		G4cout << "Now Index #" << currentIdx << G4endl;
}

void ExSingleParticleSource::SetVerbosity(int vL)
{
	verbosityLevel = vL;
	posGenerator[currentIdx]->SetVerbosity(vL);
	angGenerator[currentIdx]->SetVerbosity(vL);
	eneGenerator[currentIdx]->SetVerbosity(vL);
	timGenerator[currentIdx]->SetVerbosity(vL);

	particleList[currentParticleIdx]->timGenerator->SetVerbosity(vL);
	// G4cout << "Verbosity Set to: " << verbosityLevel << G4endl;
}

void ExSingleParticleSource::SetParticleDefinition(
    G4ParticleDefinition *aParticleDefinition)
{
	// G4cout<<"currentIdx: "<<currentIdx<<" "<<aParticleDefinition->GetParticleName()<<G4endl;
	// particle_definition[currentIdx] = aParticleDefinition;
	// particle_charge[currentIdx] = particle_definition[currentIdx]->GetPDGCharge();

	// G4cout<<"currentParticleIdx: "<<currentParticleIdx<<G4endl;
	particleList[currentParticleIdx]->fName = aParticleDefinition->GetParticleName();
	particleList[currentParticleIdx]->fDef = aParticleDefinition;
	particleList[currentParticleIdx]->fCharge = aParticleDefinition->GetPDGCharge();
}

void ExSingleParticleSource::SetParticleDefinition(G4String name)
{
	// G4cout<<"currentIdx: "<<currentIdx<<G4endl;
	// G4cout<<"currentParticleIdx: "<<currentParticleIdx<<G4endl;
	particleList[currentParticleIdx]->fName = name;
	particleList[currentParticleIdx]->UserRad = true;
}

void ExSingleParticleSource::SetParticleDefinition(int z, int a)
{
	// G4cout<<"currentIdx: "<<currentIdx<<G4endl;
	particleList[currentParticleIdx]->RadZ = z;
	particleList[currentParticleIdx]->RadA = a;
	particleList[currentParticleIdx]->fDecaySpec = new DecaySpec(particleList[currentParticleIdx]->fName, z, a);
}

void ExSingleParticleSource::GeneratePrimaryVertex(G4Event *evt)
{

	if (nParticlesToBeGen == 0)
		return;

	// If the generator contains series decay (such as Bi214), nParticlesToBeGen will > 0.
	// They have the same position but maybe not the same timestamp.
	// G4cout<<"nParticlesToBeGen: "<<nParticlesToBeGen<<G4endl;
	for (int cI = 0; cI < nParticlesToBeGen; cI++)
	{
		if (verbosityLevel >= 1)
			G4cout << "Correlated events: " << nParticlesToBeGen << " " << cI << G4endl;
		// G4cout<<"particle name: "<<particleList[cI]->fName<<G4endl;

		// Position stuff
		if (cI == 0)
			particleList[cI]->fPosition = particleList[cI]->posGenerator->GenerateOne();
		else
			particleList[cI]->fPosition = particleList[0]->fPosition;
		G4ThreeVector fPosition = particleList[cI]->fPosition;

		// Angular stuff
		particleList[cI]->fMomentumDirection = particleList[cI]->angGenerator->GenerateOne();
		G4ParticleMomentum fMomentum = particleList[cI]->fMomentumDirection;

		// Time stuff
		// An independent particle (such as Bi214 beta) obeys a uniform distribution, but a delayed alpha obeys a exp distribution.
		G4double tmp_t = particleList[cI]->timGenerator->GenerateOne(particleList[cI]->fDef) / CLHEP::second;

		if (particleList[cI]->timGenerator->bDependency == false)
		{
			if (cI == 0)
				tmp_t *= particleList[cI]->timGenerator->fSegmentSize / (1 * CLHEP::second);
		}
		//TimeStamp tmp_stamp(tmp_t);

		// Series particles. Add the time distribution to the last timestamp
		//if (cI > 0)
		//{
		//	// G4cout<<"Pre: "<<particleList[cI-1]->fTimeStamp<<G4endl;
		//	tmp_stamp.Add(particleList[cI - 1]->fTimeStamp);
		//}
		//particleList[cI]->fTimeStamp.Copy(tmp_stamp);
		//if (verbosityLevel >= 1)
		//	particleList[cI]->fTimeStamp.Print();

		// create a new vertex
		G4PrimaryVertex *vertex = new G4PrimaryVertex(particleList[cI]->fPosition, 0);
		JPSimPrimaryVertexInformation *vinfo = new JPSimPrimaryVertexInformation(particleList[cI]->RadZ, particleList[cI]->RadA);
		vertex->SetUserInformation(vinfo);

		if (particleList[cI]->UserRad)
		{
			std::vector<G4PrimaryParticle *> particle;
			int ppindex = -1;
			G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

			if (verbosityLevel >= 1)
				G4cout << "Generating Radioactive Background: " << particleList[cI]->fName << G4endl;

			G4ParticleDefinition *pardef = particleTable->FindParticle("geantino");

			double ek;
			std::vector<double> eG;
			DecayType type = particleList[cI]->fDecaySpec->GetVertex(ek, eG);
			switch (type)
			{
			case INVALID:
				return;
			case BetaPlus:
				pardef = particleTable->FindParticle("e+");
				break;
			case BetaMinus:
				pardef = particleTable->FindParticle("e-");
				break;
			case EC:
				pardef = particleTable->FindParticle("geantino");
				break;
			case Alpha:
				pardef = particleTable->FindParticle("alpha");
				break;
			default:
				break;
			}

			// G4cout<<"Generate: "<<pardef->GetParticleName()<<G4endl;

			// Energy stuff
			particleList[cI]->fEnergy = ek * MeV;
			G4double fEnergy = particleList[cI]->fEnergy;
			G4double fCharge = particleList[cI]->fCharge;
			G4ThreeVector fPolarization = particleList[cI]->fPolarization;
			// G4cout<<"fEnergy: "<<fEnergy<<G4endl;

			// create new primaries and set them to the vertex
			// No particles are generated when EC happens.
			if (type != EC)
			{
				particle.push_back(new G4PrimaryParticle(pardef));
				ppindex++;
				//fTotalTimeStamp.push_back(tmp_stamp);
				SetCurrentNewPrimary(particle[ppindex], currentIdx, pardef,
						     fPosition, fEnergy, fMomentum, fCharge, fPolarization);
				vertex->SetPrimary(particle[ppindex]);
			}

			// G4cout<<"Number of gammas: "<<eG.size()<<G4endl;
			for (size_t i = 0; i < eG.size(); i++)
			{
				particle.push_back(new G4PrimaryParticle(particleTable->FindParticle("gamma")));
				ppindex++;
				fMomentum = particleList[cI]->angGenerator->GenerateOne();
				fEnergy = eG[i] * MeV;
				fCharge = particleList[cI]->fCharge;
				G4ThreeVector fPolarization = particleList[cI]->fPolarization;
				//fTotalTimeStamp.push_back(tmp_stamp);
				SetCurrentNewPrimary(particle[ppindex], currentIdx, pardef,
						     fPosition, fEnergy, fMomentum, fCharge, fPolarization);
				vertex->SetPrimary(particle[ppindex]);
			}
		}

		else
		{
			// G4cout<<"Not a radioactive decay!"<<G4endl;
			//  Energy stuff
			particleList[cI]->fEnergy = particleList[cI]->eneGenerator->GenerateOne(particleList[cI]->fDef);
			// Time stuff
			//fTotalTimeStamp.push_back(tmp_stamp);
			// G4cout<<"Generate: "<<particleList[cI]->fDef->GetParticleName()<<" "<<tmp_stamp<<G4endl;

			// if(particleList[cI]->fDef==nullptr)
			//     particleList[cI]->fDef = G4Geantino::Geantino()
			G4PrimaryParticle *particle =
			    new G4PrimaryParticle(particleList[cI]->fDef);

			SetCurrentNewPrimary(particle, cI);
			vertex->SetPrimary(particle);
		}
		evt->AddPrimaryVertex(vertex);
	}

	// G4cout<<fTotalTimeStamp.size()<<G4endl;

	// now pass the weight to the primary vertex. CANNOT be used here!
	//  vertex->SetWeight(particle_weight);
	if (verbosityLevel > 1)
		G4cout << " Primary Vetex generated !" << G4endl;
}

void ExSingleParticleSource::SetCurrentNewPrimary(G4PrimaryParticle *particle, G4int index, G4ParticleDefinition *fDef, G4ThreeVector fPosition, G4double fEnergy, G4ParticleMomentum fMomentumDirection, G4double fCharge, G4ThreeVector fPolarization)
{
	G4double mass = fDef->GetPDGMass();
	particle->SetKineticEnergy(fEnergy);
	particle->SetMass(mass);
	particle->SetMomentumDirection(fMomentumDirection);
	particle->SetCharge(fCharge);
	particle->SetPolarization(fPolarization.x(), fPolarization.y(), fPolarization.z());
	if (verbosityLevel > 0)
	{
		G4cout << "Particle name: "
		       << fDef->GetParticleName() << G4endl;
		G4cout << "       Energy: " << fEnergy << G4endl;
		//particleList[index]->fTimeStamp.Print("");
		G4cout << "     Position: " << fPosition << G4endl;
		G4cout << "    Direction: " << fMomentumDirection
		       << G4endl;
	}
	// Set bweight equal to the multiple of all non-zero weights
	G4double fWeight = particleList[index]->eneGenerator->GetWeight() * biasRndm->GetBiasWeight();
	// pass it to primary particle
	particle->SetWeight(fWeight);
}

void ExSingleParticleSource::SetCurrentNewPrimary(G4PrimaryParticle *particle, G4int index)
{

	G4double mass = particleList[index]->fDef->GetPDGMass();
	particle->SetKineticEnergy(particleList[index]->fEnergy);
	particle->SetMass(mass);
	particle->SetMomentumDirection(particleList[index]->fMomentumDirection);
	particle->SetCharge(particleList[index]->fCharge);
	particle->SetPolarization(particleList[index]->fPolarization.x(),
				  particleList[index]->fPolarization.y(),
				  particleList[index]->fPolarization.z());
	if (verbosityLevel > 0)
	{
		G4cout << "Particle name: "
		       << particleList[index]->fDef->GetParticleName() << G4endl;
		G4cout << "       Energy: " << particleList[index]->fEnergy << G4endl;
		//particleList[index]->fTimeStamp.Print("");
		G4cout << "     Position: " << particleList[index]->fPosition << G4endl;
		G4cout << "    Direction: " << particleList[index]->fMomentumDirection
		       << G4endl;
	}
	// Set bweight equal to the multiple of all non-zero weights
	particleList[index]->fWeight = particleList[index]->eneGenerator->GetWeight() * biasRndm->GetBiasWeight();
	// pass it to primary particle
	particle->SetWeight(particleList[index]->fWeight);
}

void ExSingleParticleSource::SetRadZA(int Z, int A)
{
	RadZ = Z;
	RadA = A;
}
