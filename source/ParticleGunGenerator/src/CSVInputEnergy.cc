#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "TRandom.h"
#include "CSVInputEnergy.hh"

#include <sys/stat.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

CSVInputEnergy::CSVInputEnergy()
{
    fPGType = "Neutron Generator";
    fClassName = "CSVInputEnergy";
    fParticle = "gamma";
    fParticleEnergy = 0. * keV;
    fParticlePosition = G4ThreeVector(0., 0., 1e-8);
    fParticlePolarization = G4ThreeVector(0., 0., 0.);
    fParticleMomentumDirection = G4ThreeVector(0., 0., 1.0);
    fTypeFlag = 1001;
    fEvtFlag = 0;

    fCSVPath = "";
    fGraph = NULL;
    fRandomShoot = true;
}

CSVInputEnergy::~CSVInputEnergy()
{
    if(fGraph != NULL)
        delete fGraph;
}

bool CSVInputEnergy::IfFileExist(G4String FilePath)
{
    struct stat buffer;   
    return (stat(FilePath.data(), &buffer) == 0);
}

void CSVInputEnergy::LoadCSVFile()
{
    std::vector<std::pair<Double_t, Double_t>> EnergyMapVector;
    std::ifstream InputCSV(fCSVPath.Data());
    std::string Line;
    while (std::getline(InputCSV, Line))
    {
        std::string Number;
        std::istringstream ReadLine(Line);
        std::pair<Double_t, Double_t> EnergyVSLumi;
        std::getline(ReadLine, Number, ',');
        EnergyVSLumi.first = std::atof(Number.c_str());
        std::getline(ReadLine, Number, ',');
        EnergyVSLumi.second = std::atof(Number.c_str());
        EnergyMapVector.push_back(EnergyVSLumi);
    }
    G4cout << "Load " << EnergyMapVector.size() << " energy point" << G4endl;
    std::sort(EnergyMapVector.begin(), EnergyMapVector.end(), CSVInputEnergy::EnergySort);
    fNumOfPoint = EnergyMapVector.size();
    if(fNumOfPoint > 1000)
    {
        G4cerr << fClassName << ": Error!!! Input number of energy points large than 1000. Pleace enlaege the array size!" << G4endl;
        exit(EXIT_FAILURE);
    }
    for(Int_t iN = 0; iN < fNumOfPoint; ++iN)
    {
        fEnergy[iN]     = EnergyMapVector[iN].first;
        fLuminosity[iN] = EnergyMapVector[iN].second;
    }
    fGraph = new TGraph(fNumOfPoint, fEnergy, fLuminosity);
    fGraph->ComputeRange(fLoEnergy, fLoLumi, fHiEnergy, fHiLumi);
}

void CSVInputEnergy::Initialize(std::vector<std::string> PGParameters)
{
    Int_t NumOfPars = PGParameters.size();
    for (Int_t iP = 0; iP < NumOfPars; ++iP)
    {
        std::string::size_type IfPar = PGParameters[iP].find("-");
        if(IfPar != std::string::npos)
        {
            if (PGParameters[iP] == "-n")
            {
                if(NumOfPars > iP+1)
                {
                    fParticle = PGParameters[iP+1];
                    ++iP;
                }
            }
        }
        else
        {
            fCSVPath = PGParameters[iP];
            if (!IfFileExist(fCSVPath.Data()))
            {
                G4cerr << fClassName << ": Error!!! Input CVS file \" " << fCSVPath << " \" not exist!" << G4endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    if (fCSVPath == "") 
    {
        G4cerr << fClassName << ": Please input CVS path." << G4endl;
        exit(EXIT_FAILURE);
    }
    this->LoadCSVFile();
}

void CSVInputEnergy::Generate()
{
    Double_t energy = gRandom->Uniform(fLoEnergy, fHiEnergy);
    while(1)
    {
        if(fGraph->Eval(energy) < gRandom->Uniform(fHiLumi))
        {
            energy = gRandom->Uniform(fLoEnergy, fHiEnergy);
        }
        else
        {
            break;
        }
    }

    fParticleEnergy = energy * keV;
    if(fRandomShoot == true)
    {
        std::pair<G4ThreeVector, G4ThreeVector> PosAndMomDir = this->GenPosAndMomDir();
        fParticlePosition = PosAndMomDir.first;
        fParticleMomentumDirection = PosAndMomDir.second;
    }
    
    ++fEvtFlag;
}
