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
// MODULE:       ExGeneralParticleSource.cc
//
// Version:      2.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
// Documentation avaialable at http://reat.space.qinetiq.com/gps
//   These include:
//       User Requirement Document (URD)
//       Software Specification Documents (SSD)
//       Software User Manual (SUM): on-line version available
//       Technical Note (TN) on the physics and algorithms
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 2.0, 05/02/2004, Fan Lei, Created.
//    After changes to version 1.1 as in Geant4 v6.0
//     - Mutilple particle source definition
//     - Re-structured commands
//     - Split the task into smaller classes
//
//     - old commonds have been retained for backward compatibility, will be
//       removed in the future. 
//
// 2016, Linyan, Modified
// To include the correlated events in GeneratePrimaryVertex
// Delayed event will copy the prompt event position by default
// New position can be manually defined in macfile
///////////////////////////////////////////////////////////////////////////////
//
#include "G4Event.hh"
#include "Randomize.hh"
#include "ExGeneralParticleSource.hh"
//#include "JPSimEventInfo.hh"

ExGeneralParticleSource *ExGeneralParticleSource::fInstance = NULL;
G4bool ExGeneralParticleSource::fixTime = true; //Debug

ExGeneralParticleSource::ExGeneralParticleSource()
: multiple_vertex(false), flat_sampling(false)
{
	sourceVector.clear();
	sourceIntensity.clear();
	sourceProbability.clear();
	currentSource = new ExSingleParticleSource(theMessenger);
	sourceVector.push_back(currentSource);
	sourceIntensity.push_back(1.);
	currentSourceIdx = G4int(sourceVector.size() - 1);
	theMessenger = new ExGeneralParticleSourceMessenger(this);
	theMessenger->SetParticleGun(currentSource);
	IntensityNormalization();
	time_range=1e9;//ns
}

ExGeneralParticleSource::~ExGeneralParticleSource()
{
	delete currentSource;
	delete theMessenger;
	//G4cout<<"Debug:Hi , this is delete func"<<G4endl;
}

void ExGeneralParticleSource::AddaSource(G4double aV)
{
	//G4cout<<"Add a source"<<G4endl;
	currentSource = new ExSingleParticleSource(theMessenger);
	theMessenger->SetParticleGun(currentSource);
	sourceVector.push_back(currentSource);
	sourceIntensity.push_back(aV);
	currentSourceIdx = G4int(sourceVector.size() - 1);
	if(currentSourceIdx>0){
		currentSource->GetPosDist()->CopyFromSPSP(sourceVector[currentSourceIdx-1]->GetPosDist());
		currentSource->GetCurrentPosDist()->CopyFromSPSP(sourceVector[currentSourceIdx-1]->GetCurrentPosDist());
		//G4cout<<currentSourceIdx<<" "<<currentSource->GetPosDist()->GetConfineVolume()<<G4endl;
		currentSource->SetVerbosity(sourceVector[currentSourceIdx-1]->GetVerbosityLevel());
	}
	IntensityNormalization();
}

void ExGeneralParticleSource::IntensityNormalization()
{
	G4double total  = 0.;
	size_t i = 0 ;
	for (i = 0; i < sourceIntensity.size(); i++) 
		total += sourceIntensity[i] ;
	//
	sourceProbability.clear();
	std::vector <G4double> sourceNormalizedIntensity;
	sourceNormalizedIntensity.clear();

	sourceNormalizedIntensity.push_back(sourceIntensity[0]/total);
	sourceProbability.push_back(sourceNormalizedIntensity[0]);

	for ( i = 1 ;  i < sourceIntensity.size(); i++) {
		sourceNormalizedIntensity.push_back(sourceIntensity[i]/total);
		sourceProbability.push_back(sourceNormalizedIntensity[i] + sourceProbability[i-1]);
	}

	// set source weights here based on sampling scheme (analog/flat) and intensities
	for ( i = 0 ;  i < sourceIntensity.size(); i++) {
		if(currentSource->GetVerbosityLevel()>=1)G4cout<<"sourceIntensity.size: "<<sourceIntensity.size()<<" "<<i<<G4endl;
		if (!flat_sampling) {
			sourceVector[i]->GetBiasRndm()->SetIntensityWeight(1.);
		} else {
			sourceVector[i]->GetBiasRndm()->SetIntensityWeight(sourceNormalizedIntensity[i]*sourceIntensity.size());
		}
	}

	normalised = true;
} 

void ExGeneralParticleSource::ListSource()
{
	G4cout << " The number of particle sources is " << sourceIntensity.size() << G4endl;
	for (size_t i = 0 ; i < sourceIntensity.size(); i++)
		G4cout << "   source " << i << " intensity is " << sourceIntensity[i] << G4endl;
}

void ExGeneralParticleSource::SetCurrentSourceto(G4int aV)
{
	size_t id = size_t (aV) ;
	if ( id <= sourceIntensity.size() ) {
		currentSourceIdx = aV;
		currentSource = sourceVector[id];
		theMessenger->SetParticleGun(currentSource);
		//
	} else {
		G4cout << " source index is invalid " << G4endl;
		G4cout << "    it shall be <= " << sourceIntensity.size() << G4endl;
	}
}

void ExGeneralParticleSource::SetCurrentSourceIntensity(G4double aV)
{
	sourceIntensity[currentSourceIdx] = aV;
	normalised = false;
}

void ExGeneralParticleSource::ClearAll()
{
	currentSourceIdx = -1;
	currentSource = 0;
	sourceVector.clear();
	sourceIntensity.clear();
	sourceProbability.clear();
}

void ExGeneralParticleSource::DeleteaSource(G4int aV)
{
	size_t id = size_t (aV) ;
	if ( id <= sourceIntensity.size() ) {
		sourceVector.erase(sourceVector.begin()+aV);
		sourceIntensity.erase(sourceIntensity.begin()+aV);
		normalised = false ;
		if (currentSourceIdx == aV ) { 
			if ( sourceIntensity.size() > 0 ) { 
				currentSource = sourceVector[0];
				currentSourceIdx = 1;
			} else {
				currentSource = 0;
				currentSourceIdx = -1;
			}
		}	  		
	} else {
		G4cout << " source index is invalid " << G4endl;
		G4cout << "    it shall be <= " << sourceIntensity.size() << G4endl;
	}
} 

void ExGeneralParticleSource::GeneratePrimaryVertex(G4Event* evt)
{
	//evt->SetUserInformation(new JPSimEventInfo());
	//JPSimEventInfo* testInfo(static_cast< JPSimEventInfo*>(evt->GetUserInformation()));
	G4int nvertex, totalvertex=0;

	//G4cout<<"Enter GeneratePrimaryVertex()"<<G4endl;
	//G4cout<<"sourceIntensity.size(): "<<sourceIntensity.size()<<G4endl;

	for(size_t i=1;i<sourceIntensity.size();i++)
	{
		if(sourceVector[i]->GetVerbosityLevel()>=1){
			G4cout<<"Number of Sources: "<<sourceIntensity[i]*(GetTimeRange()*1e-9)<<" "<<sourceIntensity.size()<<" "<<i<<G4endl;
		}
		if(fixTime)
			nvertex = 1;
		else
			nvertex = CLHEP::RandPoisson::shoot(sourceIntensity[i]*(GetTimeRange()*1e-9));

		
		//G4cout<<"Debug:sourceIntensity--->"<<sourceIntensity[i]<<G4endl;
		//G4cout<<"Debug:sourceProbabily--->"<<sourceProbability[i]<<G4endl;
		for(int j=0;j<nvertex;j++)
        {
			if(sourceVector[i]->GetVerbosityLevel()>=1)
				G4cout<<"Nvertex: "<<nvertex<<" "<<j<<G4endl;
			

			sourceVector[i]->GeneratePrimaryVertex(evt);
			//G4cout<<"Debug:Hi, I create a evt"<<G4endl;
			//G4cout<<"GetTotalParticleIdx: "<<sourceVector[i]->GetTotalParticleIdx()<<G4endl;
			//for(int k=0;k<sourceVector[i]->GetTotalParticleNumber();k++)
   //         {
			//	if(fixTime)
			//		testInfo->SetEventStamp(TimeStamp(0,1000));
			//	else
			//		testInfo->SetEventStamp(sourceVector[i]->GetCurrentParticleTimeStamp(k));
			//	//G4cout<<sourceVector[i]->GetCurrentParticleTimeStamp(k)<<G4endl;
			//	totalvertex++;
			//}
			//sourceVector[i]->ClearTimeStampList();
		}
	}
}
