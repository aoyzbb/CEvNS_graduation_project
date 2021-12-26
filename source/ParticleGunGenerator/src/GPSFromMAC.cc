#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "GPSFromMAC.hh"

GPSFromMAC::GPSFromMAC()
{
    fPGType = "GPS from MAC file";
    fClassName = "GPSFromMAC";
    fParticle = "gamma";
    fParticleEnergy = 0. * MeV;
    fParticlePosition = G4ThreeVector(0., 0., 0.);
    fParticlePolarization = G4ThreeVector(0., 0., 0.);
    fParticleMomentumDirection = G4ThreeVector(0., 0., 0.);
    fTypeFlag = 10001;
    fEvtFlag = 0;

}

GPSFromMAC::~GPSFromMAC(){}

void GPSFromMAC::Initialize(std::vector<std::string> )
{}

void GPSFromMAC::Generate()
{
    ++fEvtFlag;
}