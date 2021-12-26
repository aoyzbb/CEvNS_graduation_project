
//*********************************************
//  This is Template of Cosmic Ray Generator
//                                  author:Ruiting
//

#ifndef _PGGenerator_H_
#define _PGGenerator_H_

#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

class PGGenerator
{
public:
    PGGenerator();
    ~PGGenerator();

    virtual void Initialize(std::vector<std::string>);
    virtual void Generate();

    virtual inline void SetParticle(G4String Particle){fParticle = Particle;}
    virtual inline void SetParticleEnergy(G4double ParticleEnergy){fParticleEnergy = ParticleEnergy;}
    virtual inline void SetParticlePosition(G4ThreeVector ParticlePosition){fParticlePosition = ParticlePosition;}
    virtual inline void SetParticlePolarization(G4ThreeVector ParticlePolarization){fParticlePolarization = ParticlePolarization;}
    virtual inline void SetParticleMomentumDirection(G4ThreeVector ParticleMomentumDirection){fParticleMomentumDirection = ParticleMomentumDirection;}
    virtual inline void SetTypeFlag(G4int TypeFlag){fTypeFlag = TypeFlag;}

    virtual inline G4String GetPGType(){return fPGType;}
    virtual inline G4String GetParticle(){return fParticle;}
    virtual inline G4double GetParticleEnergy(){return fParticleEnergy;}
    virtual inline G4ThreeVector GetParticlePosition(){return fParticlePosition;}
    virtual inline G4ThreeVector GetParticlePolarization(){return fParticlePolarization;}
    virtual inline G4ThreeVector GetParticleMomentumDirection(){return fParticleMomentumDirection;}
    virtual inline G4int GetTypeFlag(){return fTypeFlag;}
    virtual inline G4int GetGenerateNumber(){return fGenNumber;}

    G4double ConvertStrToDouble(std::string);
    //std::pair<G4ThreeVector, G4ThreeVector> GenPosAndMomDir();
    inline static bool mySort(double i, double j){return (i < j);}

private:
    G4String fPGType;
    G4String fClassName;
    G4String fParticle;
    G4double fParticleEnergy; //Unit: MeV
    G4ThreeVector fParticlePosition;
    G4ThreeVector fParticlePolarization;
    G4ThreeVector fParticleMomentumDirection;
    G4int fTypeFlag;
    G4int fGenNumber;
    G4int fEvtFlag;
};

#endif