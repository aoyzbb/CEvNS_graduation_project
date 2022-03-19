#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "PGGenerator.hh"

PGGenerator::PGGenerator()
{
    fBoxLength = 30;
    fBallRadius = 500;
}
PGGenerator::~PGGenerator(){}
void PGGenerator::Initialize(std::vector<std::string> PGParameters){}
void PGGenerator::Generate(){}

G4double PGGenerator::ConvertStrToDouble(std::string Str)
{
    if (!(std::isdigit((Str.c_str()[0])) || (Str.c_str()[0] == '-')))
        throw "ConvertStrToDouble: Input str is not a vaild number!";
    return std::stod(Str);
}
G4ThreeVector PGGenerator::ConvertStrToV3(std::string Str1, std::string Str2, std::string Str3)
{
    return G4ThreeVector(ConvertStrToDouble(Str1), ConvertStrToDouble(Str2), ConvertStrToDouble(Str3));
}

std::pair<G4ThreeVector, G4ThreeVector> PGGenerator::GenPosAndMomDir()
{
    G4double theta = acos(2*(G4UniformRand()-0.5));
    G4double phi = 2 * G4UniformRand() * M_PI;
    G4double my_x = fBoxLength * (G4UniformRand() - 0.5);
    G4double my_y = fBoxLength * (G4UniformRand() - 0.5);

    std::pair<G4ThreeVector, G4ThreeVector> PosAndMomDir;

    PosAndMomDir.first = G4ThreeVector(-my_x * sin(phi) + my_y * cos(theta) * cos(phi) + fBallRadius * sin(theta) * cos(phi),
                                      my_x * cos(phi) + my_y * cos(theta) * sin(phi) + fBallRadius * sin(theta) * sin(phi),
                                      -my_y * sin(theta) + fBallRadius * cos(theta)) * cm;
    PosAndMomDir.second = G4ThreeVector(-sin(theta) * cos(phi), -sin(theta) * sin(phi), -cos(theta));
    return PosAndMomDir;
}