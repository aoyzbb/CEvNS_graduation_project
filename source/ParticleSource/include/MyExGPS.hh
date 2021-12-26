//*********************************************
//  This is auto generated by G4gen 0.6
//                                  author:Qian

#ifndef _MyExGPS_H_
#define _MyExGPS_H_

#include "G4VUserPrimaryGeneratorAction.hh"
#include "PGGenerator.hh"
#include "PGGeneratorList.hh"
#include "ExGeneralParticleSource.hh"
#include "ExGeneralParticleSourceMessenger.hh"
#include "globals.hh"

class ExGeneralParticleSource;
class G4Event;

class MyExGPS : public G4VUserPrimaryGeneratorAction
{
public:
    MyExGPS();
    //MyExGPS(PGPSConfig);
    ~MyExGPS();

    virtual void GeneratePrimaries(G4Event *anEvent);

private:
    ExGeneralParticleSource *fParticleGun;
    PGGenerator* fPGGenerator;
    PGGeneratorList* fPGGeneratorList;
};

#endif
