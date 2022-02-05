//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#include "MyG4BasedAnalysis.hh"
//#include "g4root.hh"
#include "G4AnalysisManager.hh"
#include "Verbose.hh"

#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"

#include "G4ProcessType.hh"
//fNotDefined, fTransportation, fElectromagnetic, fOptical, fHadronic, fPhotolepton_hadron,
//fDecay, fGeneral, fParameterisation, fUserDefined, fParallel, fPhonon, fUCN
#include "G4DecayProcessType.hh"
//DECAY, DECAY_WithSpin, DECAY_PionMakeSpin, DECAY_Radioactive, DECAY_Unknown, DECAY_External
#include "G4HadronicProcessType.hh"
//fHadronElastic, fHadronInelastic, fCapture, fFission, fHadronAtRest, fLeptonAtRest, fChargeExchange, fRadioactiveDecay
#include "G4TransportationProcessType.hh"
//TRANSPORTATION, COUPLED_TRANSPORTATION, STEP_LIMITER, USER_SPECIAL_CUTS, NEUTRON_KILLER
#include "G4StepStatus.hh"
//fWorldBoundary, fGeomBoundary, fAtRestDoItProc, fAlongStepDoItProc, fPostStepDoItProc, fUserDefinedLimit, fExclusivelyForcedProc, fUndefined
#include "G4TrackStatus.hh"
//fAlive, fStopButAlive, fStopAndKill, fKillTrackAndSecondaries, fSuspend, fPostponeToNextEvent

#include "G4LogicalVolumeStore.hh"

#include "PGGeneratorList.hh"
#include "MyDetectorConstruction.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TH1F.h"
#pragma GCC diagnostic pop

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyG4BasedAnalysis::MyG4BasedAnalysis()
{
    SetFileName("g4output.root");
    fNumOfEvents = 0;
    fOutputLevel = 1;

    //-------
    //#ANALYSIS 1. 初始化变量
    if (fOutputLevel == 1)
    {
        engdep0 = 0;
        engdep1 = 0;
        engdep2 = 0;
        engdep3 = 0;
    }
    else
    {
        gamma_countppu = 0;
        gamma_countppd = 0;
        gamma_countpnu = 0;
        gamma_countpnd = 0;
        gamma_countnnu = 0;
        gamma_countnnd = 0;
        gamma_countnpu = 0;
        gamma_countnpd = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyG4BasedAnalysis::~MyG4BasedAnalysis()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyG4BasedAnalysis::BeginOfRunAction()
{
    if (!active)
        return;

    if (verbose)
        G4cout << "====>MyG4BasedAnalysis::BeginOfRunAction()" << G4endl;

    auto analysisManager = G4AnalysisManager::Instance();

    // Default settings
    analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

    analysisManager->SetVerboseLevel(1);
    analysisManager->OpenFile(fFileName);

    fPGGeneratorList = PGGeneratorList::GetInstance();
    fPGGenerator = fPGGeneratorList->GetGenerator();

    //-------
    //#ANALYSIS 2. 定义Ntuple结构

    // Creating 1D histograms
    //analysisManager->SetFirstHistoId(1);
    //analysisManager->CreateH1("phEng", "photon energy", 50, 0., 100); // h1 Id = 0

    // Creating 2D histograms
    //analysisManager->CreateH2("HitOnAnode", "Cherenkov photon hits on the anodes", // h2 Id = 0
    //                          50, -1000., 1000, 50, -1000., 1000.);

    // Creating ntuple
    //
    analysisManager->SetFirstNtupleId(1);
    analysisManager->CreateNtuple("particleshits", "Hits"); // ntuple Id = 1
    if (fOutputLevel == 1)
    {
        analysisManager->CreateNtupleDColumn("engpp");
        analysisManager->CreateNtupleDColumn("engpn");
        analysisManager->CreateNtupleDColumn("engnn");
        analysisManager->CreateNtupleDColumn("engnp");
        analysisManager->CreateNtupleDColumn("TruthEnergy");
        analysisManager->CreateNtupleDColumn("TruthPosX");
        analysisManager->CreateNtupleDColumn("TruthPosY");
        analysisManager->CreateNtupleDColumn("TruthPosZ");
        analysisManager->CreateNtupleDColumn("TruthMomDirX");
        analysisManager->CreateNtupleDColumn("TruthMomDirY");
        analysisManager->CreateNtupleDColumn("TruthMomDirZ");
    }
    else
    {
        analysisManager->CreateNtupleDColumn("countppu");
        analysisManager->CreateNtupleDColumn("countppd");
        analysisManager->CreateNtupleDColumn("countpnu");
        analysisManager->CreateNtupleDColumn("countpnd");
        analysisManager->CreateNtupleDColumn("countnnu");
        analysisManager->CreateNtupleDColumn("countnnd");
        analysisManager->CreateNtupleDColumn("countnpu");
        analysisManager->CreateNtupleDColumn("countnpd");
        analysisManager->CreateNtupleDColumn("TruthPosX");
        analysisManager->CreateNtupleDColumn("TruthPosY");
        analysisManager->CreateNtupleDColumn("TruthPosZ");
        analysisManager->CreateNtupleDColumn("TruthMomDirX");
        analysisManager->CreateNtupleDColumn("TruthMomDirY");
        analysisManager->CreateNtupleDColumn("TruthMomDirZ");
    }
    analysisManager->FinishNtuple();

    return;
}

void MyG4BasedAnalysis::EndOfRunAction()
{
    if (!active)
        return;

    if (verbose)
        G4cout << "====>MyG4BasedAnalysis::EndOfRunAction()" << G4endl;

    //-------
    //#ANALYSIS 6. 在Run结束的时候将ntuple保存到文件

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
     G4cout << "\n====> In total " << fNumOfEvents << " Events have been stored." << G4endl;

    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyG4BasedAnalysis::BeginOfEventAction(const G4Event *evt)
{
    if (!active)
        return;

    if (verbose > 1)
        G4cout << "====>MyG4BasedAnalysis::BeginOfEventAction()" << G4endl;

    if(fNumOfEvents%100==0)
        G4cout << "====>Generate " << fNumOfEvents << " events." << G4endl;

    //-------
    fEvent = evt;

    //#ANALYSIS 3. 初始化Event开始的参数
    if (fOutputLevel == 1)
    {
        engdep0 = 0;
        engdep1 = 0;
        engdep2 = 0;
        engdep3 = 0;
    }
    else
    {
        gamma_countppu = 0;
        gamma_countppd = 0;
        gamma_countpnu = 0;
        gamma_countpnd = 0;
        gamma_countnnu = 0;
        gamma_countnnd = 0;
        gamma_countnpu = 0;
        gamma_countnpd = 0;
    }

    return;
}

void MyG4BasedAnalysis::EndOfEventAction(const G4Event *)
{
    if (!active)
        return;

    if (verbose > 1)
        G4cout << "====>MyG4BasedAnalysis::EndOfEventAction()" << G4endl;

    //-------
    //#ANALYSIS 5. 在Event结束的时候将数据保存到ntuple

    auto analysisManager = G4AnalysisManager::Instance();
    if (fOutputLevel == 1)
    {
        analysisManager->FillNtupleDColumn(1, 0, engdep0);
        analysisManager->FillNtupleDColumn(1, 1, engdep1);
        analysisManager->FillNtupleDColumn(1, 2, engdep2);
        analysisManager->FillNtupleDColumn(1, 3, engdep3);
    }
    else
    {
        analysisManager->FillNtupleDColumn(1, 0, gamma_countppu);
        analysisManager->FillNtupleDColumn(1, 1, gamma_countppd);
        analysisManager->FillNtupleDColumn(1, 2, gamma_countpnu);
        analysisManager->FillNtupleDColumn(1, 3, gamma_countpnd);
        analysisManager->FillNtupleDColumn(1, 4, gamma_countnnu);
        analysisManager->FillNtupleDColumn(1, 5, gamma_countnnd);
        analysisManager->FillNtupleDColumn(1, 6, gamma_countnpu);
        analysisManager->FillNtupleDColumn(1, 7, gamma_countnpd);
    }
    //保存Truth信息
    fTruthEnergy = fPGGenerator->GetParticleEnergy();
    fTruthPosX = fPGGenerator->GetParticlePosition()[0];
    fTruthPosY = fPGGenerator->GetParticlePosition()[1];
    fTruthPosZ = fPGGenerator->GetParticlePosition()[2];
    fTruthMomDirX = fPGGenerator->GetParticleMomentumDirection()[0];
    fTruthMomDirY = fPGGenerator->GetParticleMomentumDirection()[1];
    fTruthMomDirZ = fPGGenerator->GetParticleMomentumDirection()[2];

    Int_t TruthOutColumn = 4;
    if (fOutputLevel != 1)
        TruthOutColumn = 8;
    analysisManager->FillNtupleDColumn(1, TruthOutColumn    , fTruthEnergy);
    analysisManager->FillNtupleDColumn(1, TruthOutColumn + 1, fTruthPosX);
    analysisManager->FillNtupleDColumn(1, TruthOutColumn + 2, fTruthPosY);
    analysisManager->FillNtupleDColumn(1, TruthOutColumn + 3, fTruthPosZ);
    analysisManager->FillNtupleDColumn(1, TruthOutColumn + 4, fTruthMomDirX);
    analysisManager->FillNtupleDColumn(1, TruthOutColumn + 5, fTruthMomDirY);
    analysisManager->FillNtupleDColumn(1, TruthOutColumn + 6, fTruthMomDirZ);

    analysisManager->AddNtupleRow(1); 
    ++fNumOfEvents;

    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ClassificationOfNewTrack MyG4BasedAnalysis::ClassifyNewTrack(const G4Track *aTrack)
{
    if (!active)
        return fUrgent;

    if (verbose > 2)
        G4cout << "====>MyG4BasedAnalysis::ClassifyNewTrack()" << G4endl;

    //-------
    //#ANALYSIS 4.1 在生成新Track的时候保存相应数据

    return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyG4BasedAnalysis::PreTrackingAction(const G4Track *aTrack)
{
    if (!active)
        return;

    if (verbose > 2)
        G4cout << "====>MyG4BasedAnalysis::PreTrackingAction()" << G4endl;

    //-------
    //#ANALYSIS 4.2 在Tracking产生的时候保存相应数据

    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyG4BasedAnalysis::PostTrackingAction(const G4Track *aTrack)
{
    if (!active)
        return;

    if (verbose > 2)
        G4cout << "====>MyG4BasedAnalysis::PostTrackingAction()" << G4endl;

    //-------
    //#ANALYSIS 4.3 在Tracking终止的时候保存相应数据

    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyG4BasedAnalysis::SteppingAction(const G4Step *aStep)
{
    if (!active)
        return;

    if (verbose > 2)
        G4cout << "====>MyG4BasedAnalysis::SteppingAction()" << G4endl;

    //---
    //1. 相关参数的获取

    //1.1 Track的相关参数
    const G4Track *aTrack = aStep->GetTrack();

    //1.2 Step的相关参数
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    //1.3 以下是拿DetectorConstruction相关参数的方法，根据具体情况修改即可

    //-------
    //#ANALYSIS 4.4 在Steppinging的时候保存相应数据

    //---
    //2. 添加一些判断，并保存对应的数据。以下为演示，且按ANALYSIS 2. Ntuple定义的结构进行保存

    G4ThreeVector postPos = postStepPoint->GetPosition();

    auto *pVolume = postStepPoint->GetTouchableHandle()->GetVolume();
    if (pVolume == NULL)
        return;

    //ntuple1: hits kernel?
    G4LogicalVolume *presentVolume = pVolume->GetLogicalVolume();
    //G4cout<<"Debug:Hi, I hit!---->"<<presentVolume->GetName()<<"<----Now I am in"<<G4endl;

    if (fOutputLevel == 1)
    {
        if (presentVolume->GetName() == "CsIvol")
        {
            G4double engdep_inside_this_step = aStep->GetTotalEnergyDeposit();
            // G4cout<<"Debug:Hi, I hit!---->"<<engdep_inside_this_step<<"<----my energy deposit"<<G4endl;

            if (postPos.x() > 0 && postPos.y() > 0)
            {
                engdep0 += engdep_inside_this_step;
            }

            if (postPos.x() > 0 && postPos.y() < 0)
            {
                engdep1 += engdep_inside_this_step;
            }

            if (postPos.x() < 0 && postPos.y() < 0)
            {
                engdep2 += engdep_inside_this_step;
            }

            if (postPos.x() < 0 && postPos.y() > 0)
            {
                engdep3 += engdep_inside_this_step;
            }
            if (postPos.x() * postPos.y() == 0)
            {
                G4cout << "Debug:error found!,x,y is" << postPos.x() << "," << postPos.y() << G4endl;
            }
        }
    }
    else
    {
        G4ParticleDefinition *particle = aTrack->GetDefinition();
        G4String name = particle->GetParticleName();
        if (postStepPoint->GetProcessDefinedStep()->GetProcessName() == "opticalphoton" && presentVolume->GetName() == "PMTvol") //要求是光子且击中PMT
        {

            if (postPos.x() > 0 && postPos.y() > 0)
            {
                if (postPos.z() > 0)
                {
                    gamma_countppu += 1;
                }
                else
                {
                    gamma_countppd += 1;
                }
            }

            if (postPos.x() > 0 && postPos.y() < 0)
            {
                if (postPos.z() > 0)
                {
                    gamma_countpnu += 1;
                }
                else
                {
                    gamma_countpnd += 1;
                }
            }

            if (postPos.x() < 0 && postPos.y() < 0)
            {
                if (postPos.z() > 0)
                {
                    gamma_countnnu += 1;
                }
                else
                {
                    gamma_countnnd += 1;
                }
            }

            if (postPos.x() < 0 && postPos.y() > 0)
            {
                if (postPos.z() > 0)
                {
                    gamma_countnpu += 1;
                }
                else
                {
                    gamma_countnpd += 1;
                }
            }
            if (postPos.x() * postPos.y() == 0)
            {
                G4cout << "Debug:error found!,x,y is" << postPos.x() << "," << postPos.y() << G4endl;
            }

            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return;
        }
    }

    return;
}
