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

    //每个闪烁体上有两个PMT, u for upside , d for downside 如:ppn:第一象限闪烁体的上PMT
    G4double gamma_countppu;   
    G4double gamma_countppd;

    G4double gamma_countpnu;
    G4double gamma_countpnd;
    
    G4double gamma_countnnu;
    G4double gamma_countnnd;
    
    G4double gamma_countnpu;
    G4double gamma_countnpd;
    
    //四个闪烁体,分布在四个象限, p for positive , n for negative; 如pp: x>0 , y>0
    G4double engdep0; //pp
    G4double engdep1; //pn
    G4double engdep2; //nn
    G4double engdep3; //np

    //Truth 信息
    G4double fTruthEnergy;
    G4double fTruthPosX;
    G4double fTruthPosY;
    G4double fTruthPosZ;
    G4double fTruthMomDirX;
    G4double fTruthMomDirY;
    G4double fTruthMomDirZ;
    
    //经过不同层时的能量分布
    G4double eng_water;
    G4double eng_ps;
    G4double eng_contemPb;
    G4double eng_lowPb;
    G4double eng_HDPE;
    G4double eng_kernel;
    G4double eng_CsI;

    //仅针对 单个 CsI 信息
    std::vector<double> eng_deposit;
    std::vector<double> pid_deposit;
    std::vector<double> x_deposit;
    std::vector<double> y_deposit;
    std::vector<int> process_type;
    std::vector<double> pid_through;
    std::vector<double> kenergy_through;

    G4double truthenergy2;
    G4bool eflag;
    

};

#endif