//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#ifndef _MyDETECTORCONSTRUCTION_H_
#define _MyDETECTORCONSTRUCTION_H_

#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"
#include "G4FieldManager.hh"

#include "MyDetectorSettings.hh"

#include <vector>

class G4GlobalMagFieldMessenger;
class MyDetectorMessenger;
class MyMagneticField;

/// Detector construction for laoding GDML geometry

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    MyDetectorConstruction(G4String GDMLFile);
    ~MyDetectorConstruction();

    virtual G4VPhysicalVolume *Construct();
    virtual void ConstructSDandField();

    void SetReadFile(const G4String &fname);
    void SetWriteFile(const G4String &fname);
    PrimaryDescription* GetPrimaryDescription() { return fSettings->GetPrimaryDescription(); }

    //#AuxXML 2. 定义与DetSetting的接口函数
    MyDetectorSettings* GetDetSettings() { return fSettings; }
    
    //#PhysTRD 3. 检查与DetSetting的接口函数
    RadiatorDescription* GetRadiatorDescription() { return fSettings->GetRadiatorDescription(); }

private:
    G4GDMLReadStructure *fReader;
    G4GDMLWriteStructure *fWriter;
    G4GDMLParser *fParser;
    G4String fReadFile, fWriteFile;
    G4bool fWritingFlag;
    
    G4GlobalMagFieldMessenger *fMagFieldMessenger;
    MyMagneticField* fMagneticField;
    G4FieldManager* fFieldMgr;

    MyDetectorSettings *fSettings;
    MyDetectorMessenger *fDetectorMessenger;
};

#endif
