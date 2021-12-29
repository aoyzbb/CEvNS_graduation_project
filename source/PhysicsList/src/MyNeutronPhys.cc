//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#include "MyNeutronPhys.hh"

//#include "NeutronHPMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"

// Processes
#include "G4CrossSectionInelastic.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPThermalScattering.hh"

#include "G4NeutronInelasticProcess.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPCapture.hh"

#include "G4NeutronInelasticCrossSection.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronFissionProcess.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPFission.hh"

#include "G4SystemOfUnits.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGNuclNuclXsc.hh"

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"

#include "G4NeutronRadCapture.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyNeutronPhys::MyNeutronPhys(const G4String &name)
    : G4VPhysicsConstructor(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyNeutronPhys::~MyNeutronPhys()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyNeutronPhys::ConstructProcess()
{
	// Model vs. Energy:
	// MeV level: HP(High Precision Particle & Neutron), LEND,
	//  Inelastic scattering models

	G4ParticleDefinition *neutron = G4Neutron::Neutron();
	G4ProcessManager *pManager = neutron->GetProcessManager();
	G4ProcessTable *processTable = G4ProcessTable::GetProcessTable();

	// delete all neutron processes if already registered
	//
	G4VProcess *process = 0;
	process = processTable->FindProcess("hadElastic", neutron);
	if (process)
		pManager->RemoveProcess(process);
	//
	process = processTable->FindProcess("neutronInelastic", neutron);
	if (process)
		pManager->RemoveProcess(process);
	//
	process = processTable->FindProcess("nCapture", neutron);
	if (process)
		pManager->RemoveProcess(process);
	//
	process = processTable->FindProcess("nFission", neutron);
	if (process)
		pManager->RemoveProcess(process);

	//-- elastic Physics
	G4HadronElasticProcess *theElasticProcess = new G4HadronElasticProcess;

	// for fast neutron
	G4HadronElastic *elastic_neutronChipsModel = new G4ChipsElasticModel();
	elastic_neutronChipsModel->SetMinEnergy(20.0 * MeV);
	theElasticProcess->RegisterMe(elastic_neutronChipsModel);
	// dataset options: G4ChipsNeutronElasticXS (dataset from CHIPS package), G4NeutronElasticXS
	theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));

	// for slow neutron
	G4ParticleHPElastic *theElasticNeutronHP = new G4ParticleHPElastic;
	theElasticNeutronHP->SetMinEnergy(4. * eV);
	theElasticNeutronHP->SetMaxEnergy(20.0 * MeV);
	theElasticProcess->RegisterMe(theElasticNeutronHP);
	theElasticProcess->AddDataSet(new G4ParticleHPElasticData);

	// for thermal neutron
	G4ParticleHPThermalScattering *thermalNeutron = new G4ParticleHPThermalScattering();
	theElasticProcess->RegisterMe(thermalNeutron);
	theElasticProcess->AddDataSet(new G4ParticleHPThermalScatteringData());

	pManager->AddDiscreteProcess(theElasticProcess);

	//-- inelastic Physics
	G4NeutronInelasticProcess *theInelasticProcess = new G4NeutronInelasticProcess("inelastic");
	theInelasticProcess->AddDataSet(new G4BGGNucleonInelasticXS(G4Neutron::Neutron()));

	G4FTFModel * theStringModel = new G4FTFModel;
	G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
	theStringModel->SetFragmentationModel( theStringDecay );
	G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
	G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

	G4TheoFSGenerator *theFTFModel1 = new G4TheoFSGenerator("FTFP");
	theFTFModel1->SetHighEnergyGenerator(theStringModel);
	theFTFModel1->SetTransport(theCascade);
	theFTFModel1->SetMinEnergy(5.0 * GeV);
	theFTFModel1->SetMaxEnergy(100.0 * TeV);
	theInelasticProcess->RegisterMe(theFTFModel1);

	G4CascadeInterface *theBERTModel1 = new G4CascadeInterface;
	theBERTModel1->SetMinEnergy(20.0 * MeV);
	theBERTModel1->SetMaxEnergy(5.0 * GeV);
	theInelasticProcess->RegisterMe(theBERTModel1);

	G4ParticleHPInelastic *theNeutronInelasticHPModel = new G4ParticleHPInelastic;
	theNeutronInelasticHPModel->SetMinEnergy(0.0 * eV);
	theNeutronInelasticHPModel->SetMaxEnergy(20.0 * MeV);
	theInelasticProcess->RegisterMe(theNeutronInelasticHPModel);
	theInelasticProcess->AddDataSet(new G4ParticleHPInelasticData);

	pManager->AddDiscreteProcess(theInelasticProcess);

	//-- n-Capture Physics
	G4HadronCaptureProcess *theCaptureProcess = new G4HadronCaptureProcess;
	theCaptureProcess->AddDataSet(new G4ParticleHPCaptureData);

	G4NeutronRadCapture *theCaptureModel = new G4NeutronRadCapture;
	//theCaptureModel->SetMinEnergy(20.0 * MeV);
	theCaptureProcess->RegisterMe(theCaptureModel);

	pManager->AddDiscreteProcess(theCaptureProcess);

	//-- nFission
	//
	G4HadronFissionProcess *theFissionProcess = new G4HadronFissionProcess();
	theFissionProcess->AddDataSet(new G4ParticleHPFissionData());
	G4ParticleHPFission *theNeutronFissionModel = new G4ParticleHPFission();
	theFissionProcess->RegisterMe(theNeutronFissionModel);

	pManager->AddDiscreteProcess(theFissionProcess);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
