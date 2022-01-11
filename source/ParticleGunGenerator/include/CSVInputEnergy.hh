
//*********************************************
//  This is Template of Cosmic Ray Generator
//                                  author:Ruiting
//

#ifndef _CSVInputEnergy_H_
#define _CSVInputEnergy_H_

#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "PGGenerator.hh"

#include "TGraph.h"

#include <vector>
#include <string>

class CSVInputEnergy: public PGGenerator
{
public:
    CSVInputEnergy();
    ~CSVInputEnergy();

    void Initialize(std::vector<std::string>);
    void GenTimeVector();
    void Generate();
    void LoadCSVFile();
    bool IfFileExist(G4String FilePath);

    inline static bool EnergySort(std::pair<Double_t, Double_t> i, std::pair<Double_t, Double_t> j){return (i.first < j.first);}

    inline void SetParticle(G4String Particle){fParticle = Particle;}
    inline void SetParticleEnergy(G4double ParticleEnergy){fParticleEnergy = ParticleEnergy;}
    inline void SetParticlePosition(G4ThreeVector ParticlePosition){fParticlePosition = ParticlePosition;}
    inline void SetParticlePolarization(G4ThreeVector ParticlePolarization){fParticlePolarization = ParticlePolarization;}
    inline void SetParticleMomentumDirection(G4ThreeVector ParticleMomentumDirection){fParticleMomentumDirection = ParticleMomentumDirection;}
    inline void SetTypeFlag(G4int TypeFlag){fTypeFlag = TypeFlag;}

    inline G4String GetPGType(){return fPGType;}
    inline G4String GetParticle(){return fParticle;}
    inline G4double GetParticleEnergy(){return fParticleEnergy;}
    inline G4ThreeVector GetParticlePosition(){return fParticlePosition;}
    inline G4ThreeVector GetParticlePolarization(){return fParticlePolarization;}
    inline G4ThreeVector GetParticleMomentumDirection(){return fParticleMomentumDirection;}
    inline G4int GetTypeFlag(){return fTypeFlag;}
    inline G4int GetGenerateNumber(){return fGenNumber;}

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

    TString fCSVPath;
    TGraph* fGraph;
    Int_t fNumOfPoint;
    Double_t fEnergy[1000];
    Double_t fLuminosity[1000];
    Double_t fLoEnergy;
    Double_t fHiEnergy;
    Double_t fLoLumi;
    Double_t fHiLumi;

    bool fRandomShoot;

};

#endif