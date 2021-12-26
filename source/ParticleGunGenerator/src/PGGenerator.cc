#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "PGGenerator.hh"

PGGenerator::PGGenerator(){}
PGGenerator::~PGGenerator(){}
void PGGenerator::Initialize(std::vector<std::string> PGParameters){}
void PGGenerator::Generate(){}

G4double PGGenerator::ConvertStrToDouble(std::string Str)
{
    if (!(std::isdigit((Str.c_str()[0])) || (Str.c_str()[0] == '-')))
        throw "ConvertStrToDouble: Input str is not a vaild number!";
    return std::stod(Str);
}

//std::pair<G4ThreeVector, G4ThreeVector> PGGenerator::GenPosAndMomDir()
//{
//    SampleSquare square =SampleSquare();
//    G4double theta = acos(2*(G4UniformRand()-0.5));
//    G4double phi = 2 * G4UniformRand() * M_PI;
//    G4double my_x = square.GetLength() * (G4UniformRand() - 0.5);
//    G4double my_y = square.GetLength() * (G4UniformRand() - 0.5);
//
//    std::pair<G4ThreeVector, G4ThreeVector> PosAndMomDir;
//
//    PosAndMomDir.first = G4ThreeVector(-my_x * sin(phi) + my_y * cos(theta) * cos(phi) + 100 * sin(theta) * cos(phi),
//                                      my_x * cos(phi) + my_y * cos(theta) * sin(phi) + 100 * sin(theta) * sin(phi),
//                                      -my_y * sin(theta) + 100 * cos(theta)) * cm;
//    PosAndMomDir.second = G4ThreeVector(-sin(theta) * cos(phi), -sin(theta) * sin(phi), -cos(theta));
//    return PosAndMomDir;
//}