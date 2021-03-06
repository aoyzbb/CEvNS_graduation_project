//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#include "MyPhysListOp.hh"

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"

#include "G4OpBoundaryProcess.hh"

#include "G4OpWLS.hh"
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "G4GenericMessenger.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(MyPhysListOp);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPhysListOp::MyPhysListOp(G4int verbose, const G4String &name)
    : G4OpticalPhysics(verbose, name),

      fMaxNumPhotons(100),
      fMaxBetaChange(10.0),
      fYieldFactor(1.),
      fExcitationRatio(0.0),
      fProfile("delta"),
      fFiniteRiseTime(false),
      fScintillationByParticleType(false),
      fScintillationTrackInfo(false),
      fInvokeSD(true),
      fCerenkovStackPhotons(true),
      fScintillationStackPhotons(true)
{
    verboseLevel = verbose;
    //fMessenger = new G4OpticalPhysicsMessenger(this);

    for (G4int i = 0; i < kNoProcess; i++)
    {
        fProcessUse.push_back(true);
        fProcessTrackSecondariesFirst.push_back(true);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPhysListOp::~MyPhysListOp()
{
    //delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyPhysListOp::PrintStatistics() const
{
    // Print all processes activation and their parameters

    for (G4int i = 0; i < kNoProcess; i++)
    {
        G4cout << "  " << G4OpticalProcessName(i) << " process:  ";
        if (!fProcessUse[i])
        {
            G4cout << "not used" << G4endl;
        }
        else
        {
            G4cout << "used" << G4endl;
            if (i == kCerenkov)
            {
                G4cout << "    Max number of photons per step: " << fMaxNumPhotons << G4endl;
                G4cout << "    Max beta change per step:       " << fMaxBetaChange << G4endl;
                if (fProcessTrackSecondariesFirst[kCerenkov])
                {
                    G4cout << "    Track secondaries first:  activated" << G4endl;
                }
                else
                {
                    G4cout << "    Track secondaries first:  inactivated" << G4endl;
                }
            }
            if (i == kScintillation)
            {
                if (fScintillationByParticleType)
                    G4cout << "    Scintillation by Particle Type:  activated " << G4endl;
                G4cout << "    Yield factor: " << fYieldFactor << G4endl;
                G4cout << "    ExcitationRatio: " << fExcitationRatio << G4endl;
                if (fProcessTrackSecondariesFirst[kScintillation])
                {
                    G4cout << "    Track secondaries first:  activated" << G4endl;
                }
                else
                {
                    G4cout << "    Track secondaries first:  inactivated" << G4endl;
                }
            }
            if (i == kWLS)
            {
                G4cout << "     WLS process time profile: " << fProfile << G4endl;
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4OpticalPhoton.hh"

void MyPhysListOp::ConstructParticle()
{
    /// Instantiate particles.

    // optical photon
    G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#define DIR_CMDS "/process/optical"
#define GUIDANCE "Commands related to the optical physics simulation engine for "

namespace UIhelpers
{
    //Helper functions to create UI commands for processes

    template <class T>
    void buildCommands(T *proc, const char *dir, const char *guidance)
    {
        //Generic function to add a "verbose" command for a process
        G4GenericMessenger *mess = new G4GenericMessenger(proc, dir, guidance);
        G4AutoDelete::Register(mess);
        G4GenericMessenger::Command &wlscmd4 = mess->DeclareMethod("verbose", &T::SetVerboseLevel,
                                                                   "Set the verbose level");
        wlscmd4.SetParameterName("ver", true);
        wlscmd4.SetDefaultValue("1");
        wlscmd4.SetStates(G4State_Idle);
    }

    void buildCommands(G4OpWLS *op)
    {
        //Build UI commands for WLS
        G4GenericMessenger *mess = new G4GenericMessenger(op, DIR_CMDS "/wls/", GUIDANCE " for WLS process.");
        G4AutoDelete::Register(mess);
        //Here, more correctly DeclareProperty should be used, but to do that, I would need public/friend
        //access of G4GenericMessenger to private members of G4Scintillation. Thus I use the approach
        //of DeclareMethod, that has a draw back: range checking does not work
        G4GenericMessenger::Command &wlscmd1 = mess->DeclareMethod("setTimeProfile", &G4OpWLS::UseTimeProfile,
                                                                   "Set the WLS time profile (delta or exponential)");
        wlscmd1.SetParameterName("profile", false);
        wlscmd1.SetCandidates("delta exponential");
        wlscmd1.SetStates(G4State_Idle);
        buildCommands(op, DIR_CMDS "/wls/", GUIDANCE " for WLS process.");
    }

    void buildCommands(G4Scintillation *ScintillationProcess)
    {
        //Build UI commands for scintillation
        G4GenericMessenger *mess = new G4GenericMessenger(ScintillationProcess, DIR_CMDS "/scintillation/", GUIDANCE " for scintillation process.");
        G4AutoDelete::Register(mess);
        G4GenericMessenger::Command &sccmd1 = mess->DeclareMethod("setFiniteRiseTime",
                                                                  &G4Scintillation::SetFiniteRiseTime,
                                                                  "Set option of a finite rise-time for G4Scintillation - If set, the G4Scintillation process expects the user to have set the constant material property FAST/SLOWSCINTILLATIONRISETIME");
        sccmd1.SetParameterName("time", false);
        sccmd1.SetStates(G4State_Idle);

        //G4GenericMessenger::Command &sccmd2 = mess->DeclareMethod("setYieldFactor",
        //                                                          &G4Scintillation::SetScintillationYieldFactor,
        //                                                          "Set scintillation yield factor");
        //sccmd2.SetParameterName("factor", false);
        //sccmd2.SetRange("factor>=0."); //LIMITATION: w/ DeclareMethod range checking does not work
        //sccmd2.SetStates(G4State_Idle);

        //G4GenericMessenger::Command &sccmd3 = mess->DeclareMethod("setExcitationRatio",
        //                                                          &G4Scintillation::SetScintillationExcitationRatio,
        //                                                          "Set scintillation excitation ratio");
        //sccmd3.SetParameterName("ratio", false);
        //sccmd3.SetRange("ratio>=0.&&ratio<=1.");//LIMITATION: w/ DeclareMethod range checking does not work
        //sccmd3.SetStates(G4State_Idle);

        G4GenericMessenger::Command &sccmd4 = mess->DeclareMethod("setByParticleType",
                                                                  &G4Scintillation::SetScintillationByParticleType,
                                                                  "Activate/Inactivate scintillation process by particle type");
        sccmd4.SetParameterName("flag", false);
        sccmd4.SetStates(G4State_Idle);

        G4GenericMessenger::Command &sccmd5 = mess->DeclareMethod("setTrackInfo",
                                                                  &G4Scintillation::SetScintillationTrackInfo,
                                                                  "Activate/Inactivate scintillation track info");
        sccmd5.SetParameterName("flag", false);
        sccmd5.SetStates(G4State_Idle);

        G4GenericMessenger::Command &sccmd6 = mess->DeclareMethod("setTrackSecondariesFirst",
                                                                  &G4Scintillation::SetTrackSecondariesFirst,
                                                                  "Set option to track secondaries before finishing their parent track");
        sccmd6.SetParameterName("flag", false);
        sccmd6.SetStates(G4State_Idle);

        buildCommands(ScintillationProcess, DIR_CMDS "/scintillation/", GUIDANCE " for scintillation process.");
    }

    void buildCommands(G4Cerenkov *CerenkovProcess)
    {
        //BUild UI commands for cerenkov
        G4GenericMessenger *mess = new G4GenericMessenger(CerenkovProcess, DIR_CMDS "/cerenkov/", GUIDANCE " for Cerenkov process.");
        G4AutoDelete::Register(mess);
        G4GenericMessenger::Command &cecmd1 = mess->DeclareMethod("setMaxPhotons",
                                                                  &G4Cerenkov::SetMaxNumPhotonsPerStep,
                                                                  "Set maximum number of photons per step");
        cecmd1.SetParameterName("max", false);
        //cecmd1.SetRange("max>=0");//LIMITATION: w/ DeclareMethod range checking does not work
        cecmd1.SetStates(G4State_Idle);

        G4GenericMessenger::Command &cecmd2 = mess->DeclareMethod("setMaxBetaChange",
                                                                  &G4Cerenkov::SetMaxBetaChangePerStep,
                                                                  "Set maximum change of beta of parent particle per step");
        cecmd2.SetParameterName("max", false);
        //cecmd2.SetRange("max>=0.");//LIMITATION: w/ DeclareMethod range checking does not work
        cecmd2.SetStates(G4State_Idle);

        G4GenericMessenger::Command &cecmd3 = mess->DeclareMethod("setTrackSecondariesFirst",
                                                                  &G4Cerenkov::SetTrackSecondariesFirst,
                                                                  "Set option to track secondaries before finishing their parent track");
        cecmd3.SetParameterName("flag", false);
        cecmd3.SetStates(G4State_Idle);

        buildCommands(CerenkovProcess, DIR_CMDS "/cerenkov/", GUIDANCE " for Cerenkov process.");
    }
} // namespace UIhelpers

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Threading.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4AutoDelete.hh"

void MyPhysListOp::ConstructProcess()
{
    // Construct optical processes.

    if (verboseLevel > 0)
        G4cout << "MyPhysListOp:: Add Optical Physics Processes" << G4endl;

    // A vector of optical processes
    std::vector<G4VProcess *> OpProcesses;

    for (G4int i = 0; i < kNoProcess; i++)
        OpProcesses.push_back(NULL);

    // Add Optical Processes

    G4OpAbsorption *OpAbsorptionProcess = new G4OpAbsorption();
    UIhelpers::buildCommands(OpAbsorptionProcess, DIR_CMDS "/absorption/", GUIDANCE " for absorption process");
    OpProcesses[kAbsorption] = OpAbsorptionProcess;

    G4OpRayleigh *OpRayleighScatteringProcess = new G4OpRayleigh();
    UIhelpers::buildCommands(OpRayleighScatteringProcess, DIR_CMDS "/rayleigh/", GUIDANCE " for Reyleigh scattering process");
    OpProcesses[kRayleigh] = OpRayleighScatteringProcess;

    G4OpMieHG *OpMieHGScatteringProcess = new G4OpMieHG();
    UIhelpers::buildCommands(OpMieHGScatteringProcess, DIR_CMDS "/mie/", GUIDANCE " for Mie cattering process");
    OpProcesses[kMieHG] = OpMieHGScatteringProcess;

    G4OpBoundaryProcess *OpBoundaryProcess = new G4OpBoundaryProcess();
    UIhelpers::buildCommands(OpBoundaryProcess, DIR_CMDS "/boundary/", GUIDANCE " for boundary process");
    OpBoundaryProcess->SetInvokeSD(fInvokeSD);
    OpProcesses[kBoundary] = OpBoundaryProcess;

    G4OpWLS *OpWLSProcess = new G4OpWLS();
    OpWLSProcess->UseTimeProfile(fProfile);
    UIhelpers::buildCommands(OpWLSProcess);
    OpProcesses[kWLS] = OpWLSProcess;

    G4ProcessManager *pManager = 0;
    pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

    if (!pManager)
    {
        std::ostringstream o;
        o << "Optical Photon without a Process Manager";
        G4Exception("MyPhysListOp::ConstructProcess()", "",
                    FatalException, o.str().c_str());
        return;
    }

    for (G4int i = kAbsorption; i <= kWLS; i++)
    {
        if (fProcessUse[i])
        {
            pManager->AddDiscreteProcess(OpProcesses[i]);
        }
    }

    G4Scintillation *ScintillationProcess = new G4Scintillation();
    //ScintillationProcess->SetScintillationYieldFactor(fYieldFactor);
    //ScintillationProcess->SetScintillationExcitationRatio(fExcitationRatio);
    ScintillationProcess->SetFiniteRiseTime(fFiniteRiseTime);
    ScintillationProcess->SetScintillationByParticleType(fScintillationByParticleType);
    ScintillationProcess->SetScintillationTrackInfo(fScintillationTrackInfo);
    ScintillationProcess->SetTrackSecondariesFirst(fProcessTrackSecondariesFirst[kScintillation]);
    ScintillationProcess->SetStackPhotons(fScintillationStackPhotons);
    G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();
    ScintillationProcess->AddSaturation(emSaturation);
    UIhelpers::buildCommands(ScintillationProcess);
    OpProcesses[kScintillation] = ScintillationProcess;

    //#PhysOpt 2. ???MyPhyListOp???????????????????????????

    // ?????????????????????alpha????????????????????????????????????
    // scintillation process for alpha:
    G4Scintillation *theScintProcessAlpha = new G4Scintillation("Scintillation");
    theScintProcessAlpha->SetTrackSecondariesFirst(true);
    //theScintProcessAlpha->SetScintillationYieldFactor(1.1);
    //theScintProcessAlpha->SetScintillationExcitationRatio(1.0);
    theScintProcessAlpha->SetVerboseLevel(verboseLevel);

    // scintillation process for heavy nuclei
    G4Scintillation *theScintProcessNuc = new G4Scintillation("Scintillation");
    theScintProcessNuc->SetTrackSecondariesFirst(true);
    //theScintProcessNuc->SetScintillationYieldFactor(0.2);
    //theScintProcessNuc->SetScintillationExcitationRatio(1.0);
    theScintProcessNuc->SetVerboseLevel(verboseLevel);

    //--------------------------------------------------
    G4Cerenkov *CerenkovProcess = new G4Cerenkov();
    CerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotons);
    CerenkovProcess->SetMaxBetaChangePerStep(fMaxBetaChange);
    CerenkovProcess->SetTrackSecondariesFirst(fProcessTrackSecondariesFirst[kCerenkov]);
    CerenkovProcess->SetStackPhotons(fCerenkovStackPhotons);
    UIhelpers::buildCommands(CerenkovProcess);
    OpProcesses[kCerenkov] = CerenkovProcess;

    auto myParticleIterator = GetParticleIterator();
    myParticleIterator->reset();

    while ((*myParticleIterator)())
    {

        G4ParticleDefinition *particle = myParticleIterator->value();
        G4String particleName = particle->GetParticleName();

        pManager = particle->GetProcessManager();
        if (!pManager)
        {
            std::ostringstream o;
            o << "Particle " << particleName << "without a Process Manager";
            G4Exception("MyPhysListOp::ConstructProcess()", "",
                        FatalException, o.str().c_str());
            return; // else coverity complains for pManager use below
        }

        if (CerenkovProcess->IsApplicable(*particle) &&
            fProcessUse[kCerenkov])
        {
            pManager->AddProcess(CerenkovProcess);
            pManager->SetProcessOrdering(CerenkovProcess, idxPostStep);
        }
        if (ScintillationProcess->IsApplicable(*particle) &&
            fProcessUse[kScintillation])
        {
            //#PhysOpt3. ???MyPhyListOp????????????????????????????????????
            if (particle->GetParticleName() == "GenericIon")
            {
                pManager->AddProcess(theScintProcessNuc); // AtRestDiscrete
                pManager->SetProcessOrderingToLast(theScintProcessNuc, idxAtRest);
                pManager->SetProcessOrderingToLast(theScintProcessNuc, idxPostStep);
            }
            else if (particle->GetParticleName() == "alpha")
            {
                pManager->AddProcess(theScintProcessAlpha);
                pManager->SetProcessOrderingToLast(theScintProcessAlpha, idxAtRest);
                pManager->SetProcessOrderingToLast(theScintProcessAlpha, idxPostStep);
            }
            else
            {
                pManager->AddProcess(ScintillationProcess);
                pManager->SetProcessOrderingToLast(ScintillationProcess, idxAtRest);
                pManager->SetProcessOrderingToLast(ScintillationProcess, idxPostStep);
            }
        }
    }

    // Add verbose
    for (G4int i = 0; i < kNoProcess; i++)
    {
        if (fProcessUse[i])
            OpProcesses[i]->SetVerboseLevel(verboseLevel);
    }

    if (verboseLevel > 1)
        PrintStatistics();
    if (verboseLevel > 0)
        G4cout << "### " << namePhysics << " physics constructed." << G4endl;
}

void MyPhysListOp::SetScintillationYieldFactor(G4double yieldFactor)
{
    /// Set the scintillation yield factor

    fYieldFactor = yieldFactor;
    //G4Scintillation::SetScintillationYieldFactor(yieldFactor);
}

void MyPhysListOp::SetScintillationExcitationRatio(G4double excitationRatio)
{
    /// Set the scintillation excitation ratio

    fExcitationRatio = excitationRatio;
    //  G4Scintillation::SetScintillationExcitationRatio(excitationRatio);
}

void MyPhysListOp::SetMaxNumPhotonsPerStep(G4int maxNumPhotons)
{
    /// Limit step to the specified maximum number of Cherenkov photons

    fMaxNumPhotons = maxNumPhotons;
    //  G4Cerenkov::SetMaxNumPhotonsPerStep(maxNumPhotons);
}

void MyPhysListOp::SetMaxBetaChangePerStep(G4double maxBetaChange)
{
    /// Limit step to the specified maximum change of beta of the parent particle

    fMaxBetaChange = maxBetaChange;
    //  G4Cerenkov::SetMaxBetaChangePerStep(maxBetaChange);
}

void MyPhysListOp::SetWLSTimeProfile(G4String profile)
{
    /// Set the WLS time profile (delta or exponential)
    fProfile = profile;
}

//void MyPhysListOp::AddScintillationSaturation(G4EmSaturation* saturation)
//{
///// Adds Birks Saturation to the G4Scintillation Process
//  G4Scintillation::AddSaturation(saturation);
//}

void MyPhysListOp::
    SetScintillationByParticleType(G4bool scintillationByParticleType)
{
    fScintillationByParticleType = scintillationByParticleType;
    //G4Scintillation::SetScintillationByParticleType(scintillationByParticleType);
}

void MyPhysListOp::
    SetScintillationTrackInfo(G4bool scintillationTrackInfo)
{
    fScintillationTrackInfo = scintillationTrackInfo;
    //G4Scintillation::SetScintillationTrackInfo(scintillationTrackInfo);
}

void MyPhysListOp::SetTrackSecondariesFirst(G4OpticalProcessIndex index,
                                            G4bool trackSecondariesFirst)
{
    if (index >= kNoProcess)
        return;
    if (fProcessTrackSecondariesFirst[index] == trackSecondariesFirst)
        return;
    fProcessTrackSecondariesFirst[index] = trackSecondariesFirst;

    //  if ( index == kCerenkov )
    //     G4Cerenkov::SetTrackSecondariesFirst(trackSecondariesFirst);
    //  if ( index == kScintillation)
    //      G4Scintillation::SetTrackSecondariesFirst(trackSecondariesFirst);
}

void MyPhysListOp::SetFiniteRiseTime(G4bool finiteRiseTime)
{
    fFiniteRiseTime = finiteRiseTime;
    //G4Scintillation::SetFiniteRiseTime(finiteRiseTime);
}

void MyPhysListOp::SetInvokeSD(G4bool invokeSD)
{
    fInvokeSD = invokeSD;
}

void MyPhysListOp::SetCerenkovStackPhotons(G4bool stackingFlag)
{
    fCerenkovStackPhotons = stackingFlag;
}

void MyPhysListOp::SetScintillationStackPhotons(G4bool stackingFlag)
{
    fScintillationStackPhotons = stackingFlag;
}

void MyPhysListOp::Configure(G4OpticalProcessIndex index, G4bool isUse)
{
    // Configure the physics constructor to use/not use a selected process.
    // This method can only be called in PreInit> phase (before execution of
    // ConstructProcess). The process is not added to particle's process manager
    // and so it cannot be re-activated later in Idle> phase with the command
    // /process/activate.

    if (index >= kNoProcess)
        return;
    if (fProcessUse[index] == isUse)
        return;
    fProcessUse[index] = isUse;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
