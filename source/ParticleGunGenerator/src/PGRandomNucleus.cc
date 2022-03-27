#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "PGGeneratorList.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Box.hh"

#include "PGRandomNucleus.hh"

#include <string>

double Spec_anti_v_mu(double *x, double *par)
{
    const double W = 52.83;     //half muon mass
    return 12. / TMath::Power(W, 4) * TMath::Power(x[0], 2) * (W - x[0]);
}

double Spec_v_e(double *x, double *par)
{
    const double W = 52.83;
    return 6. / TMath::Power(W, 4) * TMath::Power(x[0], 2) * (W - 2./3.*x[0]);
}

double Dis_Er(double* x, double *par)
{
    const double u_mass = 931.494061;                   //atomic mass in (MeV)
    double Er = x[0];
    double Ev = par[0];
    double A = par[1];

    //calculate form factor
    const double fm2MeV = 1./197.64;                    //1fm = 1/197.63 MeV^-1
    const double a = 0.7 * fm2MeV;                      //yukawa potencial range
    double ma = A * u_mass;
    double R_A = TMath::Power(A, 1./3.) * 1.2 *fm2MeV;    //in MeV^-1
    double rho_0 = 3 * A / (4 * TMath::Pi() * TMath::Power(R_A, 3.));
    double q = TMath::Sqrt(2 * ma * Er);

    double FormFactor =  4 * TMath::Pi() * rho_0 / (A * TMath::Power(q, 3.))
            * (TMath::Sin(q * R_A) - q * R_A * TMath::Cos(q * R_A))
            * 1. / (1 + TMath::Power(a*q, 2));

    return (1. - ma * Er / (2 * TMath::Power(Ev, 2))) * TMath::Power(FormFactor, 2.);
    //Form factor will cause the CosTheta_v distribution deviate from f(x) = 1/2 * (x+1), where x = cosTheta
    //As a result, the CosTheta_N distriution will deviate from g(y) = 4*y*(1-y^2)

    //return (1. - ma * Er / (2 * TMath::Power(Ev, 2)));    //ignor form factor, can be used to examine the sampling process
                                                            //By doing this, the CosTheta_v will obey f(x) and CosTheta_N will obey g(y)
}

PGRandomNucleus::PGRandomNucleus()
{
    fPGType = "Coherent elastic neutrino-nucleus scattering";
    fClassName = "PGRandomNucleus";
    fParticle = "Cs137";
    fParticleEnergy = 0. * MeV;
    fParticlePosition = G4ThreeVector(0., 0., 0.);
    fParticlePolarization = G4ThreeVector(0., 0., 0.);
    fParticleMomentumDirection = G4ThreeVector(0., 0., 0.);
    fTypeFlag = 1002;
    fEvtFlag = 0;

    fCsINum = 0;
    fCsISize = G4ThreeVector(0., 0., 0.);
    fCsIPos.clear();

    fMapNucleus[G4String("Cs137")] = std::pair<G4int, G4long>(137, 1000551370);
    fMapNucleus[G4String("I127")]  = std::pair<G4int, G4long>(127, 1000531270);

    fRandomNucleus = new RandomNucleus();
}

PGRandomNucleus::~PGRandomNucleus()
{
    delete fRandomNucleus;
}

void PGRandomNucleus::Initialize(std::vector<std::string> PGParameters)
{
    G4int NumOfPars = PGParameters.size();
    for (G4int iP = 0; iP < NumOfPars; ++iP)
    {
        std::string::size_type IfPar = PGParameters[iP].find("-");
        if(IfPar != std::string::npos)
        {
            if (PGParameters[iP] == "-n")
            {
                if(NumOfPars > iP+1)
                {
                    fParticle = PGParameters[iP+1];
                    fRandomNucleus->SetA(fMapNucleus[fParticle].first);
                    ++iP;
                }
            }
        }
    }
    this->GetCsIVolAndPos();
}
void PGRandomNucleus::GetCsIVolAndPos()
{
    G4LogicalVolumeStore* LogVolStore = G4LogicalVolumeStore::GetInstance();
    G4LogicalVolume* SDLogVol = LogVolStore->GetVolume("CsIvol");

    //Get SD Size
    G4Box* SDSolVol = dynamic_cast<G4Box*>(SDLogVol->GetSolid());
    fCsISize[0] = SDSolVol->GetXHalfLength();
    fCsISize[1] = SDSolVol->GetYHalfLength();
    fCsISize[2] = SDSolVol->GetZHalfLength();

    fCsINum = UInt_t(SDLogVol->GetNoDaughters());
    for (UInt_t iP = 0; iP < fCsINum; ++iP)
    {
        G4VPhysicalVolume* SDPhyVol = SDLogVol->GetDaughter(iP);
        G4ThreeVector VSDPosition = SDPhyVol->GetTranslation();
        fCsIPos.push_back(VSDPosition);
    }
}

G4ThreeVector PGRandomNucleus::NucleusPos()
{
    G4ThreeVector nulcleusPos;
    nulcleusPos[0] = gRandom->Uniform(-fCsISize[0], fCsISize[0]);
    nulcleusPos[1] = gRandom->Uniform(-fCsISize[1], fCsISize[1]);
    nulcleusPos[2] = gRandom->Uniform(-fCsISize[2], fCsISize[2]);

    G4int whichCsI = gRandom->Integer(fCsINum);
    nulcleusPos = nulcleusPos + fCsIPos[whichCsI];
    return nulcleusPos;
}

G4ThreeVector PGRandomNucleus::NucleusMomDir(G4double theta)
{
    G4double phi = gRandom->Uniform(2.*TMath::Pi());
    return G4ThreeVector(TMath::Sin(theta) * TMath::Cos(phi), TMath::Sin(theta) * TMath::Sin(phi), TMath::Cos(theta));
}

void PGRandomNucleus::Generate()
{
    NucleusThetaE aNucleus = fRandomNucleus->GetRndmNucleusThetaE();
    fParticleEnergy = aNucleus.E;
    fParticlePosition = this->NucleusPos();
    fParticleMomentumDirection = this->NucleusMomDir(aNucleus.Theta);
    ++fEvtFlag;
}

RandomNucleus::RandomNucleus()
{
    f_Spec_anti_v_mu = new TF1("f_Spec_anti_v_mu", Spec_anti_v_mu, 0., W);
    f_Spec_anti_v_mu -> SetParameter(0, W);
    //f_Spec_anti_v_mu -> SetNpx(f_Npx);

    f_Spec_v_e = new TF1("f_Spec_v_e", Spec_v_e, 0., W);
    f_Spec_v_e -> SetParameter(0, W);
    //f_Spec_v_e -> SetNpx(f_Npx);

    f_Dis_Er = new TF1("f_Dis_Er", Dis_Er, 0., 100., 2);
    //f_Dis_Er -> SetNpx(f_Npx);
}

RandomNucleus::RandomNucleus(int A)
{
    f_A = A;                              //set nuclear mass number 

    //gRandom = new TRandomMixMax();
    //f_seed = time(NULL);  
    //gRandom -> SetSeed(f_seed);         //by default, use the time() as the seed

    f_Spec_anti_v_mu = new TF1("f_Spec_anti_v_mu", Spec_anti_v_mu, 0., W);
    f_Spec_anti_v_mu -> SetParameter(0, W);
    //f_Spec_anti_v_mu -> SetNpx(f_Npx);

    f_Spec_v_e = new TF1("f_Spec_v_e", Spec_v_e, 0., W);
    f_Spec_v_e -> SetParameter(0, W);
    //f_Spec_v_e -> SetNpx(f_Npx);

    f_Dis_Er = new TF1("f_Dis_Er", Dis_Er, 0., 100., 2);
    //f_Dis_Er -> SetNpx(f_Npx);
}

RandomNucleus::~RandomNucleus()
{
    delete gRandom;
    gRandom = NULL;

    delete f_Spec_anti_v_mu;
    f_Spec_anti_v_mu = NULL;
    
    delete f_Spec_v_e;
    f_Spec_v_e = NULL;

    delete f_Dis_Er;
    f_Dis_Er = NULL;
}

void RandomNucleus::SetSeed(long seed)
{
    f_seed = seed;
    gRandom -> SetSeed(f_seed);
}

long RandomNucleus::GetSeed()
{
    return f_seed;
}

void RandomNucleus::SetNpx(int Npx)         //warning: calling this function will significantly slow the sampling speed,
{                                           //Even when setting Npx to 100, same as the default value, the speed will be about 10 times slower. The reason is unkonwn       
    f_Npx = Npx;                            //since our functions are not rather smoth, the defualt setting is good enough
    f_Spec_anti_v_mu -> SetNpx(f_Npx);      //So by default we don't call this function
    f_Spec_v_e -> SetNpx(f_Npx);
    f_Dis_Er -> SetNpx(f_Npx);
}

int RandomNucleus::RandomvType()
{
    double rndm = gRandom -> Rndm();
    if(rndm <= 1./3.)
        return v_mu;
    else if(rndm > 1./3. && rndm <= 2./3.)
        return v_anti_mu;
    else
        return v_e;
}

NucleusThetaE RandomNucleus::GetRndmNucleusThetaE()
{
    int f_vType = RandomvType();
    double Ev = RandomEv(f_vType);
    double Er = RandomEr(Ev);
    double CosTheta_v = 1 - Er * f_A * u_mass / TMath::Power(Ev, 2.);   //With aproximation Er << Ev
    double CosTheta_N = TMath::Sqrt((1 - CosTheta_v) / 2.);             //With aproximation Er << Ev
    double Theta_N = TMath::ACos(CosTheta_N);

    NucleusThetaE theta_e;
    theta_e.E = Er;
    theta_e.Theta = Theta_N;
    return theta_e;
}

std::vector<NucleusThetaE> RandomNucleus::GetRndmNucleusThetaE(int N)
{
    std::vector<NucleusThetaE> N_Theta_Es(N);
    for(int n = 0; n < N; n++)
    {
        N_Theta_Es[n] = GetRndmNucleusThetaE();
    }
    return N_Theta_Es;
}

double RandomNucleus::RandomEv(int v_type)
{
    if(v_type == v_mu)
        return 29.79;
    else if(v_type == v_anti_mu)
    {
        return f_Spec_anti_v_mu -> GetRandom(0., W, static_cast<TRandom*>(gRandom));
    }
    else if(v_type == v_e)
    {
        return f_Spec_v_e -> GetRandom(0., W, static_cast<TRandom*>(gRandom));
    }
    else
    {
        return -1.;
    }
}

double RandomNucleus::RandomEr(double Ev)
{
    double Er_max = 2 * TMath::Power(Ev, 2.) / (f_A * u_mass);
    f_Dis_Er -> SetParameters(Ev, f_A);
    f_Dis_Er -> SetRange(0., Er_max);
    return f_Dis_Er -> GetRandom(0., Er_max, static_cast<TRandom*>(gRandom));
}