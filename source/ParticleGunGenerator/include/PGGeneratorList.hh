//*********************************************
//  This is List of all the Cosmic Ray
//                                  author:Ruiting
//

#ifndef _PGGeneratorList_H_
#define _PGGeneratorList_H_

#include "G4String.hh"
#include "PGGenerator.hh"

#include <map>

class PGGenerator;
class GammaRayBursts;
class PGGeneratorList;

struct PGPSConfig
{
    //Config for ParticleGun
    bool PGEnable;
    G4String ParticleGunType;
    G4double GenTime;
    G4int GenEvents;
    G4int GenValidEvents;
    bool OnlyValid;
    std::vector<std::string> ParticleGunParameters;
    G4String AdditionalMACCommand;
    //Simple Particle Gun
    G4String ParticleName;
    G4double ParticleEnergy;
    G4ThreeVector ParticlePosition;
    G4ThreeVector ParticlePolarization;
    G4ThreeVector ParticleMomentumDirection;
    //Config for GPS
    bool PSEnable;
    G4String PSSignal;
    //Config for ExGPX
    bool ExGPSEnable;
    G4String ExGPSMacFile;
};


class PGGeneratorList
{
public:
    PGGeneratorList();
    ~PGGeneratorList();

    static PGGeneratorList *GetInstance()
    {
        if (PGGeneratorList::fInstance == NULL)
            PGGeneratorList::fInstance = new PGGeneratorList();
        return PGGeneratorList::fInstance;
    }

    inline void SetPGPSConfig(PGPSConfig ConfigPS){fConfigPS = ConfigPS;}

    PGGenerator* GetGenerator();
    PGGenerator* GetGenerator(G4String);
    inline PGPSConfig GetPGPSConfig() const {return fConfigPS;}

    inline void SetValidEvent(G4bool Valid){fValid = Valid;}
    inline G4bool IfValidEvent(){return fValid;}

private:
    static PGGeneratorList* fInstance;
    PGPSConfig fConfigPS;
    G4bool fValid;
    std::map<G4String, PGGenerator*> fMapGeneratorPtr;

};

#endif