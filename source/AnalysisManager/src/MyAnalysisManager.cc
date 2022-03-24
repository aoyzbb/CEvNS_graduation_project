//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//
//
// The purpose of this class is to store the simulated data to a given nutple/root
// There are two ways:
// 1. use g4root, G4AnalysisManager, the geant4 ntuple manager to store the data.
// 2. use root, MyAnalysisManager->SimEvent or tree/branch, the root file to store the data.
// either way can work, choose one to fullfill your requirements.

#include "MyAnalysisManager.hh"
#include "MyRootBasedAnalysis.hh"

#include "MyDetectorConstruction.hh"
#include "MyG4BasedAnalysis.hh"
#include "MyRootBasedAnalysis.hh"
#include "Verbose.hh"

#include "G4EventManager.hh"
#include "G4TrackingManager.hh"

#include "PGGeneratorList.hh"

#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyAnalysisManager *MyAnalysisManager::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyAnalysisManager::MyAnalysisManager()
{
    // default using G4AnalysisManager
    fMyG4BasedAnalysis = new MyG4BasedAnalysis();
    fMyRootBasedAnalysis = new MyRootBasedAnalysis();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyAnalysisManager::~MyAnalysisManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyAnalysisManager::SetOutputInfo(OutputControl ControlOutput) 
{
    fControlOutput = ControlOutput;

    fMyG4BasedAnalysis->Activated();
    fMyG4BasedAnalysis->SetFileName(fControlOutput.OutputFile);
    fMyG4BasedAnalysis->SetOutputLevel(fControlOutput.OutputLevel);
    //fMyRootBasedAnalysis->Activated();
    //fMyRootBasedAnalysis->SetFileName(fControlOutput.OutputFile);
}

void MyAnalysisManager::BeginOfRunAction()
{
    if (verbose)
        G4cout << "====>MyAnalysisManager::BeginOfRunAction()" << G4endl;

    fMyG4BasedAnalysis->BeginOfRunAction();
    fMyRootBasedAnalysis->BeginOfRunAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyAnalysisManager::EndOfRunAction()
{
    if (verbose)
        G4cout << "====>MyAnalysisManager::EndOfRunAction()" << G4endl;

    fMyG4BasedAnalysis->EndOfRunAction();
    fMyRootBasedAnalysis->EndOfRunAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyAnalysisManager::BeginOfEventAction(const G4Event *evt)
{
    if (verbose)
        G4cout << "====>MyAnalysisManager::BeginOfEventAction()" << G4endl;
    if(fControlOutput.IfTrackVerbose)
    {
        G4int EventID = evt->GetEventID();
        if(std::find(fControlOutput.VerboseEvents.begin(), fControlOutput.VerboseEvents.end(), EventID) != fControlOutput.VerboseEvents.end())
        {
            G4EventManager *eventManager = G4EventManager::GetEventManager();
            eventManager->SetVerboseLevel(1);
            G4TrackingManager *trackManager = eventManager->GetTrackingManager();
            trackManager->SetVerboseLevel(1);
        }
    }

    fMyG4BasedAnalysis->BeginOfEventAction(evt);
    fMyRootBasedAnalysis->BeginOfEventAction(evt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyAnalysisManager::EndOfEventAction(const G4Event *evt)
{
    if (verbose)
        G4cout << "====>MyAnalysisManager::EndOfEventAction()" << G4endl;

    if(fControlOutput.IfTrackVerbose)
    {
        G4EventManager *eventManager = G4EventManager::GetEventManager();
        if(eventManager->GetVerboseLevel() != 0)
        {
            eventManager->SetVerboseLevel(0);
            G4TrackingManager *trackManager = eventManager->GetTrackingManager();
            trackManager->SetVerboseLevel(0);
        }
    }
    fMyG4BasedAnalysis->EndOfEventAction(evt);
    fMyRootBasedAnalysis->EndOfEventAction(evt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ClassificationOfNewTrack MyAnalysisManager::ClassifyNewTrack(const G4Track *aTrack)
{
    if (fMyG4BasedAnalysis->IsActivated())
        return fMyG4BasedAnalysis->ClassifyNewTrack(aTrack);
    else
        return fMyRootBasedAnalysis->ClassifyNewTrack(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyAnalysisManager::PreTrackingAction(const G4Track *aTrack)
{
    fMyG4BasedAnalysis->PreTrackingAction(aTrack);
    fMyRootBasedAnalysis->PreTrackingAction(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyAnalysisManager::PostTrackingAction(const G4Track *aTrack)
{
    fMyG4BasedAnalysis->PostTrackingAction(aTrack);
    fMyRootBasedAnalysis->PostTrackingAction(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyAnalysisManager::SteppingAction(const G4Step *aStep)
{
    fMyG4BasedAnalysis->SteppingAction(aStep);
    fMyRootBasedAnalysis->SteppingAction(aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyAnalysisManager::SetG4BasedFileName(G4String fname) { fMyG4BasedAnalysis->SetFileName(fname); }
void MyAnalysisManager::G4BasedActivated() { fMyG4BasedAnalysis->Activated(); }
void MyAnalysisManager::G4BasedDeactivated() { fMyG4BasedAnalysis->Deactivated(); }

void MyAnalysisManager::SetRootBasedFileName(G4String fname) { fMyRootBasedAnalysis->SetFileName(fname); }
void MyAnalysisManager::RootBasedActivated() { fMyRootBasedAnalysis->Activated(); }
void MyAnalysisManager::RootBasedDeactivated() { fMyRootBasedAnalysis->Deactivated(); }

void MyAnalysisManager::SetDetectorConstruction(MyDetectorConstruction *det)
{
    fDetector = det;
    fMyG4BasedAnalysis->SetDetectorConstruction(det);
    fMyRootBasedAnalysis->SetDetectorConstruction(det);
}

void MyAnalysisManager::SetRunAction(MyRunAction *run)
{
    fRunAction = run;
    fMyG4BasedAnalysis->SetRunAction(run);
    fMyRootBasedAnalysis->SetRunAction(run);
}

void MyAnalysisManager::SetEvtAction(MyEvtAction *evt)
{
    fEvtAction = evt;
    fMyG4BasedAnalysis->SetEvtAction(evt);
    fMyRootBasedAnalysis->SetEvtAction(evt);
}
