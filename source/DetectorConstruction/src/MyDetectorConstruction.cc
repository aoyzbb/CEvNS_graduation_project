//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#include "Verbose.hh"
#include "MyDetectorConstruction.hh"
#include "MySDDetector.hh"
#include "G4SDManager.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "MyDetectorReader.hh"
#include "MyDetectorMessenger.hh"
#include "MyDetectorSettings.hh"
#include "MyMagneticField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorConstruction::MyDetectorConstruction()
    : G4VUserDetectorConstruction(),
      fReader(0), fWriter(0), fParser(0), fSettings(0)
{
    if (verbose)
        G4cout << "====>MyDetectorConstruction::MyDetectorConstruction()" << G4endl;

    fReadFile = "./gdml/main.gdml";
    fWriteFile = "./gdml/main_writer.gdml";
    fWritingFlag = 0;

    fDetectorMessenger = new MyDetectorMessenger(this);

    fReader = new MyDetectorReader;
    fParser = new G4GDMLParser(fReader, fWriter);

    fSettings = new MyDetectorSettings;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorConstruction::MyDetectorConstruction(G4String GDMLFile)
    : G4VUserDetectorConstruction(),
      fReader(0), fWriter(0), fParser(0), fSettings(0)
{
    if (verbose)
        G4cout << "====>MyDetectorConstruction::MyDetectorConstruction()" << G4endl;

    fReadFile = GDMLFile;
    fWriteFile = "./gdml/main_writer.gdml";
    fWritingFlag = 0;

    fDetectorMessenger = new MyDetectorMessenger(this);

    fReader = new MyDetectorReader;
    fParser = new G4GDMLParser(fReader, fWriter);

    fSettings = new MyDetectorSettings;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorConstruction::~MyDetectorConstruction()
{
    delete fDetectorMessenger;
    delete fReader;
    delete fParser;
    delete fSettings;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    if (verbose)
        G4cout << "====>MyDetectorConstruction::Construct()" << G4endl;

    fParser->Read(fReadFile, false);
    fSettings->ApplyAuxValue(fParser);

    //Get volume by it's name...
    //G4LogicalVolumeStore::GetInstance()->GetVolume("Volume Name");

    //#Verb 2. ??????Material?????????
    //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    return fParser->GetWorldVolume();
}

void MyDetectorConstruction::ConstructSDandField()
{
    if (verbose)
        G4cout << "====>MyDetectorConstruction::ConstructSDandField()" << G4endl;

    //------------------------------------------------
    // ??????SD
    //------------------------------------------------

    /*
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    MySDRich *aRICH = new MySDRich("RICH");
    SDman->AddNewDetector(aRICH);

    //------------------------------------------------
    // Sensitive detectors GDML ??????
    //------------------------------------------------
    //  ???RICH.gdml???????????????????????????
    //        <volume name="ReadoutBoxVol">                     ??????????????????
    //      	    <materialref ref="G4_Cu" />                 ???????????????????????????
    //      		<solidref ref="ReadoutBox" />               ?????????????????????????????????
    //      		<auxiliary auxtype="SensDet" auxvalue="RICH"/>    SensDet????????????????????????????????????RICH?????????????????????????????????SDMananger???????????????
    //	      </volume>
    //------------------------------------------------
    // ???????????????RICH???logicVolume????????????????????????????????????????????????????????????PhysicsVolume
    //
    const G4GDMLAuxMapType *auxmap = fParser->GetAuxMap();
    for (G4GDMLAuxMapType::const_iterator iter = auxmap->begin(); iter != auxmap->end(); iter++)
    {
        for (G4GDMLAuxListType::const_iterator vit = (*iter).second.begin(); vit != (*iter).second.end(); vit++)
        {
            if ((*vit).type == "SensDet")
            {
                G4cout << "Attaching sensitive detector " << (*vit).value
                       << " to volume " << ((*iter).first)->GetName()
                       << " " << ((*iter).first) << G4endl << G4endl;

                G4VSensitiveDetector *mydet = SDman->FindSensitiveDetector((*vit).value);
                if (mydet)
                {
                    G4LogicalVolume *myvol = (*iter).first;
                    myvol->SetSensitiveDetector(mydet);
                }
            }
        }
    }
    */

    //------------------------------------------------
    // Fields

    //1. ???????????? ------------------------------------------------
    //#MagField 1. ??????????????????
    // ??????XML??????????????????????????????????????????????????????B=0???????????????mac???????????????????????????????????????
    //       /globalField/setValue 1.0 0 0 tesla

    if (fSettings->GetUserMagFieldFlag() == 0)
    {
        G4ThreeVector fieldValue;
        fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
        fMagFieldMessenger->SetVerboseLevel(1);

        // Register the field messenger for deleting
        G4AutoDelete::Register(fMagFieldMessenger);
    }

    //2. ???????????? ------------------------------------------------
    //#MagField 2. ??????????????????
    // ???????????????????????? MagTubeVol ?????????????????????????????????????????????????????? XML ??????????????????

    if (fSettings->GetUserMagFieldFlag() != 0)
    {
        fMagneticField = new MyMagneticField();
        fFieldMgr = new G4FieldManager();
        fFieldMgr->SetDetectorField(fMagneticField);
        fFieldMgr->CreateChordFinder(fMagneticField);
        G4bool forceToAllDaughters = true;
        auto *fMagneticLogical = G4LogicalVolumeStore::GetInstance()->GetVolume("MagTubeVol");
        fMagneticLogical->SetFieldManager(fFieldMgr, forceToAllDaughters);

        // Register the field and its manager for deleting
        G4AutoDelete::Register(fMagneticField);
        G4AutoDelete::Register(fFieldMgr);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyDetectorConstruction::SetReadFile(const G4String &fname)
{
    if (verbose)
        G4cout << "====>MyDetectorConstruction::SetReadFile()" << G4endl;

    fReadFile = fname;
    fWritingFlag = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyDetectorConstruction::SetWriteFile(const G4String &fname)
{
    if (verbose)
        G4cout << "====>MyDetectorConstruction::SetWriteFile()" << G4endl;

    fWriteFile = fname;
    fWritingFlag = 1;
}
