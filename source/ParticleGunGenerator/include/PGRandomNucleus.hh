
//*********************************************
//  This is Template of Cosmic Ray Generator
//                                  author:Ruiting
//

#ifndef _PGRandomNucleus_H_
#define _PGRandomNucleus_H_

#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "PGGenerator.hh"

#include "TMath.h"
#include "TRandom.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <map>

enum vType
{
    v_mu,
    v_anti_mu,
    v_e,
};

typedef struct NucleusThetaE
{
    double Theta;   //in rad, polar angle respect to the incoming direction of neutrino
    double E;       //in MeV
}NucleusThetaE;


class RandomNucleus
{
    private:
        //TRandomMixMax *RandomGen;           //random number generator
        int f_A = 137;                        //mass number of nucleus
        const double u_mass = 931.494061;   //atomic mass in (MeV)
        long f_seed = 1;                    //random seed
        const double W = 52.83;             //MeV, half mass of muon
        int f_Npx = 1000;                   //number of intergral bins, default Npx = 1000

        TF1 *f_Spec_anti_v_mu;              //spectrum of anti_v_mu
        TF1 *f_Spec_v_e;                    //spectrum of v_e
        TF1 *f_Dis_Er;                      //Distribution of Er

    public:
        RandomNucleus();
        RandomNucleus(int A);
        ~RandomNucleus();
         
        void SetSeed(long seed);            //Set the seed manually
        long GetSeed();

        void SetNpx(int Npx);
        inline void SetA(int A) {f_A = A;}
        
        int RandomvType();
        double RandomEv(int v_type);
        double RandomEr(double Ev);

        NucleusThetaE GetRndmNucleusThetaE();
        std::vector<NucleusThetaE> GetRndmNucleusThetaE(int N);
};

class PGRandomNucleus: public PGGenerator
{
public:
    PGRandomNucleus();
    ~PGRandomNucleus();

    void Initialize(std::vector<std::string>);
    void Generate();

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

    inline std::pair<G4String, std::pair<G4int, G4int> > GetIon() {return std::pair<G4String, std::pair<G4int, G4int> >(fParticle, fMapNucleus[fParticle]);}

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

    G4ThreeVector NucleusPos();
    G4ThreeVector NucleusMomDir(G4double);
    void GetCsIVolAndPos();

    std::map<G4String, std::pair<G4int, G4int> > fMapNucleus;
    G4int fCsINum;
    G4ThreeVector fCsISize;
    std::vector<G4ThreeVector> fCsIPos;
    RandomNucleus* fRandomNucleus;
};
#endif