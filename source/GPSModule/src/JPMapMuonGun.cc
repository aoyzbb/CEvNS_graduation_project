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
#include "JPMapMuonGun.hh"
//#include "JPSimEventInfo.hh"
#include "JPSimPrimaryVertexInformation.hh"
#include "TMath.h"
#include "G4SystemOfUnits.hh"
#include <numeric>

int JPMapMuonGun::nBeta = 0;

JPMapMuonGun::JPMapMuonGun()
{
}

JPMapMuonGun::~JPMapMuonGun()
{
}

void JPMapMuonGun::GeneratePrimaryVertex(G4Event* evt)
{
    //G4cout<<evt->GetEventID()<<G4endl;
	G4double er = G4UniformRand();

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	if(er<=0.5)
		SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));
	else
		SetParticleDefinition(particleTable->FindParticle(particleName="mu-"));

    int betaSeg = nBeta;
    int cosalphaSeg = evt->GetEventID()/100%10;
    int costhetaSeg = evt->GetEventID()/10%10;
    int phiSeg  = evt->GetEventID()/1%10;
    
    double cosalphaStart = -1+0.1*cosalphaSeg;
    double betaStart = 2*TMath::Pi()/10*betaSeg;

    double costhetaStart = 1-0.1*costhetaSeg;
    double phiStart = 2*TMath::Pi()/10*phiSeg;

    G4ThreeVector muonDir, posInSci, muonPos;
    
    muonDir = G4ThreeVector(0, 0, -1);
    muonDir.setTheta(acos(cosalphaStart+0.1*G4UniformRand()));
    muonDir.setPhi(betaStart+2*TMath::Pi()/10*G4UniformRand());
    
    posInSci = G4ThreeVector(0.645*m, 0, 0);
    posInSci.setTheta(acos(costhetaStart-0.1*G4UniformRand()));
    posInSci.setPhi(phiStart+2*TMath::Pi()/10*G4UniformRand());
    // Rotate 
    G4ThreeVector axis = muonDir.cross(G4ThreeVector(0,0,-1));
    posInSci.rotate(-G4ThreeVector(0,0,-1).angle(muonDir), axis);

    muonPos = posInSci-2.*m*muonDir;
    /*
    do
    {
        double costheta = 2*G4UniformRand()-1;
        //double costheta = -0.5;
        double theta = acos(costheta);
        double phi = 2*TMath::Pi()*G4UniformRand();
        double cosalpha = G4UniformRand();
        double beta = 2*TMath::Pi()*G4UniformRand();
        posInSci = G4ThreeVector(0.645*m, 0, 0);
        posInSci.setTheta(theta);
        posInSci.setPhi(phi);
        muonDir = G4ThreeVector(1,0,0);
        muonDir.setTheta(TMath::Pi()-acos(cosalpha));
        muonDir.setPhi(beta);
        G4ThreeVector deltaPos = -muonDir;
        muonPos = posInSci + 3.*m*deltaPos;
    }while (muonDir*posInSci>0);
    */

	Double_t Energy = 300*GeV;
	SetParticleMomentumDirection(muonDir);
    SetParticleEnergy(Energy);
	SetParticlePosition(muonPos);

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
