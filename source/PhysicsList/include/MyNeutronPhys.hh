//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#ifndef MyNeutronPhys_h
#define MyNeutronPhys_h 1

#include "G4VPhysicsConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MyNeutronPhys : public G4VPhysicsConstructor
{
public:
    MyNeutronPhys(const G4String &name = "neutron");
    ~MyNeutronPhys();

    virtual void ConstructParticle() {};
    virtual void ConstructProcess();

private:
};

#endif
