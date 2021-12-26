//*********************************************
//  This is auto generated by G4gen 0.6
//                                  author:Qian

#ifndef MyStepAction_h
#define MyStepAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include "MyRunAction.hh"
#include "MyEvtAction.hh"
#include "MyDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MyStepAction : public G4UserSteppingAction
{
public:
    MyStepAction(MyRunAction *, MyEvtAction *, MyDetectorConstruction *);
    virtual ~MyStepAction();

    virtual void UserSteppingAction(const G4Step *);

private:
    MyRunAction *fRun;
    MyEvtAction *fEvt;
    MyDetectorConstruction *fDetector;
};

#endif
