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
#ifndef JPCosmicMuonGun_H
#define JPCosmicMuonGun_H 1

#include "globals.hh"
#include <vector>

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "TF2.h"
#include "TH3D.h"
#include "TFile.h"

class SeaLevelMuon;

class JPCosmicMuonGun : public G4ParticleGun
{
  //
public:

	JPCosmicMuonGun(); 
	~JPCosmicMuonGun(); 

	void GeneratePrimaryVertex(G4Event*); 
	static TString muonFile;


private:
	//TF2* f2;
	//TF1* f1;
	TFile* f;
	TH3D* h3;
	Double_t R[5];
	TH3D* h[5];
    G4double L;
    G4double W;
    G4double H;
    G4double CenterHeight;
	//SeaLevelMuon* sfptr;


private:
	//static Double_t costhetastarf(Double_t costheta);
	//static Double_t Estar(Double_t E, Double_t costhetastar);
	//static Double_t distribution(Double_t* x,Double_t* par);
	//G4ThreeVector normVectors[5];

  
};

#endif
