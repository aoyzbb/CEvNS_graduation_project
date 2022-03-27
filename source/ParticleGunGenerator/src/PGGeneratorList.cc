#include "GPSFromMAC.hh"
#include "SimpleParticleGun.hh"
#include "PGNeutron.hh"
#include "CSVInputEnergy.hh"
#include "PGRandomNucleus.hh"

#include "PGGeneratorList.hh"

PGGeneratorList *PGGeneratorList::fInstance = NULL;

PGGeneratorList::PGGeneratorList():fValid(false)
{
    fMapGeneratorPtr[G4String("GPS")]  = new GPSFromMAC();
    fMapGeneratorPtr[G4String("Simple")]  = new SimpleParticleGun();
    fMapGeneratorPtr[G4String("PGN")]  = new PGNeutron();
    fMapGeneratorPtr[G4String("CSV")]  = new CSVInputEnergy();
    fMapGeneratorPtr[G4String("Nucleus")] = new PGRandomNucleus();
}

PGGeneratorList::~PGGeneratorList(){
    for(const auto& GenPtr:fMapGeneratorPtr)
        delete GenPtr.second;
}

PGGenerator* PGGeneratorList::GetGenerator()
{
    if (fMapGeneratorPtr.find(fConfigPS.ParticleGunType) != fMapGeneratorPtr.end())
        return fMapGeneratorPtr[fConfigPS.ParticleGunType];
    else
    {
        G4cerr << "Error!!! ParticleGun Type \" " << fConfigPS.ParticleGunType << " \" NOT FOUND!!! Please check your PGType config." << G4endl;
        G4cerr << "Avaiable Type:" << G4endl;
        for(const auto& GenPtr:fMapGeneratorPtr)
            G4cerr << "  " << GenPtr.first << ": " << (GenPtr.second)->GetPGType() << G4endl;
        
        exit(EXIT_FAILURE);
    }
}

PGGenerator* PGGeneratorList::GetGenerator(G4String PGType)
{
    if (fMapGeneratorPtr.find(PGType) != fMapGeneratorPtr.end())
        return fMapGeneratorPtr[PGType];
    else
    {
        G4cerr << "Error!!! ParticleGun Type \" " << PGType << " \" NOT FOUND!!! Please check your PGType config." << G4endl;
        G4cerr << "Avaiable Type:" << G4endl;
        for(const auto& GenPtr:fMapGeneratorPtr)
            G4cerr << "  " << GenPtr.first << ": " << (GenPtr.second)->GetPGType() << G4endl;
        
        exit(EXIT_FAILURE);
    }
}