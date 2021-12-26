#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "PGNeutron.hh"

PGNeutron::PGNeutron()
{
    fPGType = "Neutron Generator";
    fClassName = "PGNeutron";
    fParticle = "neutron";
    fParticleEnergy = 0. * MeV;
    fParticlePosition = G4ThreeVector(0., 0., -750.);
    fParticlePolarization = G4ThreeVector(0., 0., 0.);
    fParticleMomentumDirection = G4ThreeVector(0., 0., 1.0);
    fTypeFlag = 1001;
    fEvtFlag = 0;

    norm = 0.0379;
}

PGNeutron::~PGNeutron(){}

void PGNeutron::Initialize(std::vector<std::string> PGParameters)
{
    if (PGParameters.size() > 0)
    {
        if(PGParameters.size() != 1)
        {
            G4cout << fClassName << ": Warning!!! Length of Input PGParameters and number of Initialize parameters of PGNeutron not matched!!!" << G4endl;
            G4cout << fClassName << ": Use default Initialize parameters" << G4endl;
        }
        else
        {
            try
            {
                norm = ConvertStrToDouble(PGParameters[0]);
            }
            catch (const char *msg)
            {
                G4cerr << msg << G4endl;
                G4cerr << fClassName << ": Please check your PGParameter config" << G4endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}

void PGNeutron::Generate()
{
    double randomnum = G4UniformRand();
    int loop = 0;
    while (randomnum <= 0.01)
    {
        randomnum = G4UniformRand();
        loop += 1;
        if (loop > 100000)
        {
            G4cout << "Warning:Can't find the energy!" << G4endl;
            break;
        }
    }
    double norm = 0.0379; //power law ,  rank = -1.6 , 5/3 A *0.01 ** (-0.6) = 1 (0.01: 10 keV)

    double Eng = pow((1 - randomnum) / norm, -5 / 3);

    fParticleEnergy = Eng * MeV;
    ++fEvtFlag;
}
