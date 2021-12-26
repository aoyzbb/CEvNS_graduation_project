//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#include "MyExGPS.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

#include "Verbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyExGPS::MyExGPS()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(0)
{
    if (verbose)
        G4cout << "====>MyExGPS::MyExGPS()" << G4endl;

    fParticleGun  = ExGeneralParticleSource::GetInstance();
    fPGGeneratorList = PGGeneratorList::GetInstance();
    fPGGenerator = fPGGeneratorList->GetGenerator();

    fPGGenerator->SetParticle(fParticleGun->GetParticleDefinition()->GetParticleName());


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyExGPS::~MyExGPS()
{
    //delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyExGPS::GeneratePrimaries(G4Event *anEvent)
{
    if (verbose)
        G4cout << "====>MyExGPS::GeneratePrimaries()" << G4endl;

    fParticleGun->GeneratePrimaryVertex(anEvent);
    fPGGenerator->SetParticleEnergy(fParticleGun->GetParticleEnergy());
    fPGGenerator->SetParticlePosition(fParticleGun->GetParticlePosition());
    fPGGenerator->SetParticlePolarization(fParticleGun->GetParticlePolarization());
    fPGGenerator->SetParticleMomentumDirection(fParticleGun->GetParticleMomentumDirection());
}