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
//
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "JPUniformMuonGun.hh"
#include "JPSimPrimaryVertexInformation.hh"
#include <numeric>
#include "ModifiedGaisserFormula.h"
#include "G4SystemOfUnits.hh"


JPUniformMuonGun::JPUniformMuonGun()
{
}

JPUniformMuonGun::~JPUniformMuonGun()
{
}

void JPUniformMuonGun::GeneratePrimaryVertex(G4Event* evt)
{
	G4double er = G4UniformRand();

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	if(er<=0.5)
		SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));
	else
		SetParticleDefinition(particleTable->FindParticle(particleName="mu-"));
	
	Double_t costheta = G4UniformRand();
	Double_t theta = acos(costheta);
	Double_t Energy = 300*GeV;
	Double_t phi = 2*TMath::Pi()*G4UniformRand();
    G4ThreeVector vertexPos(3*m, 0, 0);
    vertexPos.setTheta(theta);
    vertexPos.setPhi(phi);
	
	G4ThreeVector v(0,0,-1);
	theta = TMath::Pi()-theta;
	phi = TMath::Pi()+phi;
	v.setRThetaPhi(1, theta, phi);
	SetParticleMomentumDirection(v);
    SetParticleEnergy(Energy);
	SetParticlePosition(vertexPos);

	if(!particle_definition) 
	{
		G4ExceptionDescription ED;
		ED << "Particle definition is not defined." << G4endl;
		ED << "G4ParticleGun::SetParticleDefinition() has to be invoked beforehand." << G4endl;
		G4Exception("G4ParticleGun::GeneratePrimaryVertex()","Event0109",
					FatalException, ED);
		return;
	}

	// create a new vertex
	G4PrimaryVertex* vertex = 
		new G4PrimaryVertex(particle_position,particle_time);

	// create new primaries and set them to the vertex
	G4double mass =  particle_definition->GetPDGMass();
	for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ ){
		G4PrimaryParticle* particle =
			new G4PrimaryParticle(particle_definition);
		particle->SetKineticEnergy( particle_energy );
		particle->SetMass( mass );
		particle->SetMomentumDirection( particle_momentum_direction );
		particle->SetCharge( particle_charge );
		particle->SetPolarization(particle_polarization.x(),
					particle_polarization.y(),
					particle_polarization.z());
		vertex->SetPrimary( particle );
	}

	//evt->SetUserInformation(new JPSimEventInfo());
	//JPSimEventInfo* testInfo(static_cast< JPSimEventInfo*>(evt->GetUserInformation()));
	//testInfo->SetEventStamp(TimeStamp(0, 1000));
	JPSimPrimaryVertexInformation* vinfo = new JPSimPrimaryVertexInformation;
	vertex->SetUserInformation(vinfo);

	evt->AddPrimaryVertex( vertex );

}
