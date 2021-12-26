#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "PGGeneratorList.hh"

#include "SimpleParticleGun.hh"

SimpleParticleGun::SimpleParticleGun()
{
    fPGType = "Simple Particle Gun";
    fClassName = "SimpleParticleGun";
    fParticle = "gamma";
    fParticleEnergy = 0. * MeV;
    fParticlePosition = G4ThreeVector(0., 0., 0.);
    fParticlePolarization = G4ThreeVector(0., 0., 0.);
    fParticleMomentumDirection = G4ThreeVector(0., 0., 0.);
    fTypeFlag = 10001;
    fEvtFlag = 0;

}

SimpleParticleGun::~SimpleParticleGun(){}

void SimpleParticleGun::Initialize(std::vector<std::string> PGParameters)
{
    PGGeneratorList* CRList = PGGeneratorList::GetInstance();
    PGPSConfig Config =  CRList->GetPGPSConfig();
    fParticle = Config.ParticleName;
    fParticleEnergy = Config.ParticleEnergy * MeV;
    fParticlePosition = Config.ParticlePosition * mm;
    fParticlePolarization = Config.ParticlePolarization;
    fParticleMomentumDirection = Config.ParticleMomentumDirection;
}

void SimpleParticleGun::Generate()
{
    ++fEvtFlag;
}