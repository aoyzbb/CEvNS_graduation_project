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
#include "JPCosmicMuonGun.hh"
//#include "JPSimEventInfo.hh"
#include "JPSimPrimaryVertexInformation.hh"
#include <numeric>
#include "ModifiedGaisserFormula.h"
#include "G4SystemOfUnits.hh"

TString JPCosmicMuonGun::muonFile = "MuonUnderground.root";


JPCosmicMuonGun::JPCosmicMuonGun()
{
	G4cout<<"Muon spectrum file: "<<muonFile<<G4endl;
	f = new TFile(muonFile);
    L = 4*m;
    W = 4*m;
    H = 6*m;
    CenterHeight = 1.75*m;

	h3 = (TH3D*)f->Get("MuonDist");

	Int_t nbinsX = h3->GetNbinsX();	// Ek
	Int_t nbinsY = h3->GetNbinsY();	// costheta
    Int_t nbinsZ = h3->GetNbinsZ();	// phi

    G4cout<<nbinsX<<"\t"<<nbinsY<<"\t"<<nbinsZ<<G4endl;
	
	// R[0]: alpha = 0, beta = 0
	for(Int_t i=0; i<5; i++)
	{
		R[i] = 0;
		h[i] = (TH3D*)h3->Clone(TString::Format("h%d",i));
		h[i]->Reset();
	}
	
	for(Int_t k=1; k<=nbinsX; k++)
	{
		for(Int_t i=1; i<=nbinsY; i++)
		{
			for(Int_t j=1; j<=nbinsZ; j++)
			{
				Double_t costheta = ((TAxis*)h3->GetYaxis())->GetBinCenter(i);
				Double_t tantheta = sqrt(1-costheta*costheta)/costheta;
				Double_t phi = ((TAxis*)h3->GetZaxis())->GetBinCenter(j) + 18*degree;
				if(phi>TMath::Pi()) phi = phi-TMath::TwoPi();
				Double_t content = h3->GetBinContent(k,i,j);
				R[0] += content;
				h[0]->SetBinContent(k,i,j,content);
				if(phi>-TMath::Pi()/2 && phi<TMath::Pi()/2)
				{
					Double_t delta = content*tantheta*TMath::Cos(phi);
					R[1] += delta;
					h[1]->SetBinContent(k,i,j, delta); 
				}
				else
				{
					Double_t delta = -content*tantheta*TMath::Cos(phi);
					R[2] += delta;
					h[2]->SetBinContent(k,i,j, delta); 

				}
				if(phi>0 && phi<TMath::Pi())
				{
					Double_t delta = content*tantheta*TMath::Sin(phi);
					R[3] += delta;
					h[3]->SetBinContent(k,i,j, delta); 
				}
				else
				{
					Double_t delta = -content*tantheta*TMath::Sin(phi);
					R[4] += delta;
					h[4]->SetBinContent(k,i,j, delta); 
				}
			}

		}
	}
    

    R[0] *= W*L;
    R[1] *= H*L;
    R[2] *= H*L;
    R[3] *= H*W;
    R[4] *= H*W;

	Double_t RTotal = 0;
	for(Int_t i=0; i<5; i++)
		RTotal += R[i];

	for(Int_t i=0; i<5; i++)
	{
		G4cout<<i<<": "<<R[i]/(m*m)/h3->GetEntries()<<"m^2";
		R[i] /= RTotal;
        G4cout<<"\t"<<R[i]<<G4endl;
	}
    G4cout<<"Total area: "<<RTotal/(m*m)/h3->GetEntries()<<"m^2"<<G4endl;
    
}

JPCosmicMuonGun::~JPCosmicMuonGun()
{
	f->Close();
}

void JPCosmicMuonGun::GeneratePrimaryVertex(G4Event* evt)
{
	G4double er = G4UniformRand();

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	if(er<=0.5)
	{
		SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));
	}
	else
	{
		SetParticleDefinition(particleTable->FindParticle(particleName="mu-"));
	}
	//SetParticleDefinition(particleTable->FindParticle(particleName="geantino"));
	
	Double_t costheta = 1;
	Double_t theta = 0;
	Double_t Energy = 100;
	Double_t phi = 0;
	G4double x = 0.*m;
	G4double y = 0.*m;
	G4double z = 0.*m;
	
    
	G4double r1 = G4UniformRand();
	G4double flag1 = 0;
	size_t idx1 = 0;
	for(size_t i=0; i<5; i++)
	{
		flag1 += R[i];
		if(r1<flag1)
		{
			idx1 = i;
			break;
		}
	}
	h[idx1]->GetRandom3(Energy, costheta, phi);
	Energy = Energy*GeV;

	//G4double boxLength = 3.*m;
	// Sprinkle muon on this surface
	switch(idx1)
	{
	case 0:
	{
		x = (2*G4UniformRand()-1)*W/2;
		y = (2*G4UniformRand()-1)*L/2;
		z = H/2+CenterHeight;
		break;
	}
	case 1:
	{
		y = (2*G4UniformRand()-1)*L/2;
		z = (2*G4UniformRand()-1)*H/2+CenterHeight;
		x = W/2;
		break;
	}
	case 2:
	{
		y = (2*G4UniformRand()-1)*L/2;
		z = (2*G4UniformRand()-1)*H/2+CenterHeight;
		x = -W/2;
		break;
	}
	case 3:
	{
		x = (2*G4UniformRand()-1)*W/2;
		z = (2*G4UniformRand()-1)*H/2+CenterHeight;
		y = L/2;
		break;
	}
	case 4:
	{
		x = (2*G4UniformRand()-1)*W/2;
		z = (2*G4UniformRand()-1)*H/2+CenterHeight;
		y = -L/2;
		break;
	}
	default:
		break;

	}
	
	G4ThreeVector v(0,0,-1);
	theta = TMath::Pi()-TMath::ACos(costheta);
	phi = TMath::Pi()+phi;
	v.setRThetaPhi(1, theta, phi);
	SetParticleMomentumDirection(v);
    SetParticleEnergy(Energy);
	SetParticlePosition(G4ThreeVector(x,y,z));

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
