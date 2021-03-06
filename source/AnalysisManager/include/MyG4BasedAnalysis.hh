//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#ifndef MyG4basedAnalysis_h
#define MyG4basedAnalysis_h 1

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include <iomanip>
#include <vector>

#include "PGGeneratorList.hh"

class TH1F;
class MyDetectorConstruction;
class MyEvtAction;
class MyRunAction;

class MyG4BasedAnalysis
{
public:
    MyG4BasedAnalysis();
    ~MyG4BasedAnalysis();

    void BeginOfEventAction(const G4Event *evt);
    void EndOfEventAction(const G4Event *evt);

    void BeginOfRunAction();
    void EndOfRunAction();

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *aTrack);
    void SteppingAction(const G4Step *aStep);
    void PreTrackingAction(const G4Track *aTrack);
    void PostTrackingAction(const G4Track *aTrack);

    inline void SetFileName(G4String fname) { fFileName = fname; }
    inline void SetOutputLevel(G4int val) { fOutputLevel = val; }
    inline void Activated() { active = true; }
    inline void Deactivated() { active = false; }
    inline bool IsActivated() { return active; }

    inline void SetDetectorConstruction(MyDetectorConstruction *det) { fDetector = det; }
    inline void SetRunAction(MyRunAction *run) { fRunAction = run; }
    inline void SetEvtAction(MyEvtAction *evt) { fEvtAction = evt; }

private:
    bool active = false;
    G4String fFileName;
    MyDetectorConstruction *fDetector;
    MyRunAction *fRunAction;
    MyEvtAction *fEvtAction;
    G4int fNumOfEvents;
    PGGenerator* fPGGenerator;
    PGGeneratorList* fPGGeneratorList;

    const G4Event *fEvent;

    //#ANALYSIS 0. 定义用户变量
private:

    G4int fOutputLevel;

    G4double EnergyPP;
    G4double EnergyPN;
    G4double EnergyNN;
    G4double EnergyNP;

    //Truth 信息
    G4double fTruthEnergy;
    G4double fTruthPosX;
    G4double fTruthPosY;
    G4double fTruthPosZ;
    G4double fTruthMomDirX;
    G4double fTruthMomDirY;
    G4double fTruthMomDirZ;
    


};

#endif