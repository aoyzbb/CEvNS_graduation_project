//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#include "Verbose.hh"
#include "MyActionInitialization.hh"
#include "MyParticleGun.hh"
#include "MyPrimaryGeneratorGPS.hh"
#include "MyRunAction.hh"
#include "MyEvtAction.hh"
#include "MyStepAction.hh"
#include "MyStackAction.hh"
#include "MyTrackAction.hh"
#include "MyExGPS.hh"
#include "PGGeneratorList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyActionInitialization::MyActionInitialization(MyDetectorConstruction *detConstruction)
    : G4VUserActionInitialization(),
      fDetectorConstruction(detConstruction)
{
    if (verbose)
        G4cout << "====>MyActionInitialization::MyActionInitialization()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyActionInitialization::~MyActionInitialization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyActionInitialization::BuildForMaster() const
{
    // No specific action for Master
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyActionInitialization::Build() const
{
    if (verbose)
        G4cout << "====>MyActionInitialization::Build()" << G4endl;

    PGGeneratorList* CRList = PGGeneratorList::GetInstance();
    if ((CRList->GetPGPSConfig()).PSEnable)
    {
        G4cout << "---> Set to GPS mode." << G4endl;
        SetUserAction(new MyPrimaryGeneratorGPS());
    }
    else if((CRList->GetPGPSConfig()).PGEnable)
    {
        G4cout << "---> Set to Single Particle mode." << G4endl;
        SetUserAction(new MyParticleGun(fDetectorConstruction));
    }
    else if((CRList->GetPGPSConfig()).ExGPSEnable)
    {
        G4cout << "---> Set to ExGPS mode." << G4endl;
        SetUserAction(new MyExGPS());
    }

    MyRunAction *runAction = new MyRunAction(fDetectorConstruction);
    SetUserAction(runAction);

    MyEvtAction *evtAction = new MyEvtAction(runAction);
    SetUserAction(evtAction);

    SetUserAction(new MyStepAction(runAction, evtAction, fDetectorConstruction));

    SetUserAction(new MyStackAction());

    SetUserAction(new MyTrackAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
