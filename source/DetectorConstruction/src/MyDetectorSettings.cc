//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#include "MyDetectorSettings.hh"

#include "G4NistManager.hh"
#include "G4GDMLAuxStructType.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TROOT.h"
#include "TColor.h"
#pragma GCC diagnostic pop

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorSettings::MyDetectorSettings()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorSettings::~MyDetectorSettings()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyDetectorSettings::ApplyAuxValue(G4GDMLParser *fParser)
{
    //apply auxvalue to each logical volume
    const G4LogicalVolumeStore *lvs = G4LogicalVolumeStore::GetInstance();
    std::vector<G4LogicalVolume *>::const_iterator lvciter;
    for (lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++)
    {
        G4GDMLAuxListType auxInfo = fParser->GetVolumeAuxiliaryInformation(*lvciter);

        if (auxInfo.size() > 0)
            G4cout << "Auxiliary Information is found for Logical Volume :  "
                   << (*lvciter)->GetName() << G4endl;

        ApplyAuxValue(&auxInfo, (*lvciter));
    }

    //apply global auxvalues
    G4cout << "\nApply global auxiliary settings(if any):" << G4endl;
    ApplyAuxValue(fParser->GetAuxList());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyDetectorSettings::ApplyAuxValue(const G4GDMLAuxListType *auxInfoList, G4LogicalVolume *vol)
{
    for (std::vector<G4GDMLAuxStructType>::const_iterator iaux = auxInfoList->begin();
         iaux != auxInfoList->end(); iaux++)
    {
        G4String type = iaux->type;
        G4String value = iaux->value;
        G4String unit = iaux->unit;

        G4cout << " | " << type << " : " << value << " " << unit << G4endl;

        //#AuxXML 3. DetSetting?????????type????????????
        if (type == "UserMag")
            UserMagField = 1;

        if (type == "setColor")
            setColor(vol, value, unit);
        if (type == "setAlpha")
            setAlpha(vol, value);

        if (type == "setStepLimit")
            setStepLimit(vol, value, unit);

        //#PhysTRD 4. ??????DetSetting??????????????????
        if (type == "SetFoilThickness")
            SetFoilThickness(vol, value, unit);
        if (type == "SetGasThickness")
            SetGasThickness(vol, value, unit);
        if (type == "SetFoilNumber")
            SetFoilNumber(vol, value);
        if (type == "SetFoilMaterial")
            SetFoilMaterial(vol, value);
        if (type == "SetGasMaterial")
            SetGasMaterial(vol, value);

        //#RegCtrl 4. ??????DetSeting??????????????????
        if (type == "DefReg")
            DefineRegion(vol, value);

        if (iaux->auxList)
        {
            G4cout << "List: | " << type << " : " << value << " " << unit << G4endl;
            if (type == "PrimaryGen")
                SetPrimaryGen(iaux->auxList, value);
            else if (type == "RegionCut")
                SetRegionCut(iaux->auxList, value);
            else
                ApplyAuxValue(iaux->auxList, vol);
        }
    }
    return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//#AuxXML 4. DetSetting?????????????????????
//????????????
void MyDetectorSettings::setColor(G4LogicalVolume *vol, G4String value, G4String unit)
{
    if (vol == NULL)
        return;
    TColor *color;
    int cid = kRed;

    //TCanvas c1; [View]->[Color]??????root??????????????????
    if (value == "kRed" || value == "red")
        cid = kRed;
    if (value == "kPink" || value == "pink")
        cid = kPink;
    if (value == "kMagenta" || value == "magenta")
        cid = kMagenta;
    if (value == "kViolet" || value == "violet")
        cid = kViolet;
    if (value == "kBlue" || value == "blue")
        cid = kBlue;
    if (value == "kAzure" || value == "azure")
        cid = kAzure;
    if (value == "kCyan" || value == "cyan")
        cid = kCyan;
    if (value == "kTeal" || value == "teal")
        cid = kTeal;
    if (value == "kGreen" || value == "green")
        cid = kGreen;
    if (value == "kSpring" || value == "spring")
        cid = kSpring;
    if (value == "kYellow" || value == "yellow")
        cid = kYellow;
    if (value == "kOrange" || value == "orange")
        cid = kOrange;
    if (value == "kGray" || value == "gray")
        cid = kGray;
    if (value == "kWhite" || value == "white")
        cid = kWhite;
    if (value == "kBlack" || value == "black")
        cid = kBlack;

    color = gROOT->GetColor(cid + atoi(unit));

    if (color == NULL)
        return;

    G4VisAttributes *attrPtr = new G4VisAttributes(G4Color(color->GetRed(), color->GetGreen(), color->GetBlue(), 0.8));
    vol->SetVisAttributes(attrPtr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyDetectorSettings::setAlpha(G4LogicalVolume *vol, G4String value)
{
    if (vol == NULL)
        return;

    const G4VisAttributes *attrPtr = vol->GetVisAttributes();
    if (attrPtr == NULL)
        return;

    G4Colour color = attrPtr->GetColor();
    vol->SetVisAttributes(new G4VisAttributes(G4Color(color.GetRed(), color.GetGreen(), color.GetBlue(), atof(value))));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyDetectorSettings::setStepLimit(G4LogicalVolume *vol, G4String value, G4String unit)
{
    if (vol == NULL)
        return;

    G4UserLimits *fStepLimit = new G4UserLimits(atof(value) * G4UnitDefinition::GetValueOf(unit));
    vol->SetUserLimits(fStepLimit);
}

void MyDetectorSettings::SetPrimaryGen(const G4GDMLAuxListType *auxInfoList, G4String gunType)
{
    G4cout<<"---->gunType="<<gunType<<G4endl;
    if (gunType == "ParticleGun")
    {
        if (fPrimaryDescription == 0)
            fPrimaryDescription = new PrimaryDescription();

        for (std::vector<G4GDMLAuxStructType>::const_iterator iaux = auxInfoList->begin();
             iaux != auxInfoList->end(); iaux++)
        {
            G4String type = iaux->type;
            G4String value = iaux->value;
            G4String unit = iaux->unit;
            G4double numb = atof(value) * ((unit=="") ? 1 : G4UnitDefinition::GetValueOf(unit));

            G4cout << " |--> " << type << " : " << value << " " << unit << G4endl;

            if (type == "particle")
                fPrimaryDescription->particle = value;
            if (type == "positionX")
                fPrimaryDescription->particlePos.setX(numb);
            if (type == "positionY")
                fPrimaryDescription->particlePos.setY(numb);
            if (type == "positionZ")
                fPrimaryDescription->particlePos.setZ(numb);
            if (type == "momentumX")
                fPrimaryDescription->particleMom.setX(numb);
            if (type == "momentumY")
                fPrimaryDescription->particleMom.setY(numb);
            if (type == "momentumZ")
                fPrimaryDescription->particleMom.setZ(numb);
            if (type == "polarizationX")
                fPrimaryDescription->particlePol.setX(numb);
            if (type == "polarizationY")
                fPrimaryDescription->particlePol.setY(numb);
            if (type == "polarizationZ")
                fPrimaryDescription->particlePol.setZ(numb);
        }
    }
    else 
    {
        if (fPrimaryDescription == 0)
            fPrimaryDescription = new PrimaryDescription();

        fPrimaryDescription->PGorGPS = 1.;
        G4cout<<"Setting to GPS"<<G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//#PhysTRD 5. ??????DetSetting????????????????????????
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyDetectorSettings::SetFoilThickness(G4LogicalVolume *vol, G4String value, G4String unit)
{
    if (fRadiatorDescription == 0)
        fRadiatorDescription = new RadiatorDescription();

    fRadiatorDescription->fLogicalVolume = vol;
    fRadiatorDescription->fFoilThickness = atof(value) * G4UnitDefinition::GetValueOf(unit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyDetectorSettings::SetGasThickness(G4LogicalVolume *vol, G4String value, G4String unit)
{
    if (fRadiatorDescription == 0)
        fRadiatorDescription = new RadiatorDescription();

    fRadiatorDescription->fLogicalVolume = vol;
    fRadiatorDescription->fGasThickness = atof(value) * G4UnitDefinition::GetValueOf(unit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyDetectorSettings::SetFoilNumber(G4LogicalVolume *vol, G4String value)
{
    if (fRadiatorDescription == 0)
        fRadiatorDescription = new RadiatorDescription();

    fRadiatorDescription->fLogicalVolume = vol;
    fRadiatorDescription->fFoilNumber = atoi(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyDetectorSettings::SetFoilMaterial(G4LogicalVolume *vol, G4String value)
{
    if (fRadiatorDescription == 0)
        fRadiatorDescription = new RadiatorDescription();

    G4Material *material = G4NistManager::Instance()->FindOrBuildMaterial(value);
    if (material == 0)
        return;

    fRadiatorDescription->fLogicalVolume = vol;
    fRadiatorDescription->fFoilMaterial = material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MyDetectorSettings::SetGasMaterial(G4LogicalVolume *vol, G4String value)
{
    if (fRadiatorDescription == 0)
        fRadiatorDescription = new RadiatorDescription();

    G4Material *material = G4NistManager::Instance()->FindOrBuildMaterial(value);
    if (material == 0)
        return;

    fRadiatorDescription->fLogicalVolume = vol;
    fRadiatorDescription->fGasMaterial = material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//#RegCtrl 5. ??????DetSetting????????????????????????
void MyDetectorSettings::DefineRegion(G4LogicalVolume *vol, G4String value)
{
    G4Region *fRegion = new G4Region(value);
    vol->SetRegion(fRegion);
    fRegion->AddRootLogicalVolume(vol);
}

// RegionCut????????????????????????stopping range????????????
// ???G4?????????????????????track?????????material??????cut range????????????????????????????????????????????????????????????????????????
// ?????????????????????????????????tracker??????????????????cut?????????0.1mm??????0.01mm???????????????????????????????????????????????????????????????
// ?????????????????????????????????EMcalorimeter??????????????????shower?????????????????????????????????/????????????cut?????????0.1mm??????1mm
void MyDetectorSettings::SetRegionCut(const G4GDMLAuxListType *auxInfoList, G4String regName)
{
    G4Region *region = G4RegionStore::GetInstance()->GetRegion(regName);

    G4ProductionCuts *cuts = new G4ProductionCuts;
    for (std::vector<G4GDMLAuxStructType>::const_iterator iaux = auxInfoList->begin();
         iaux != auxInfoList->end(); iaux++)
    {
        G4String type = iaux->type;
        G4String value = iaux->value;
        G4String unit = iaux->unit;
        G4double numb = atof(value) * G4UnitDefinition::GetValueOf(unit);

        G4cout << " |--> " << type << " : " << value << " " << unit << G4endl;

        if (type == "gammaCut")
            cuts->SetProductionCut(numb, G4ProductionCuts::GetIndex("gamma"));
        if (type == "e-Cut")
            cuts->SetProductionCut(numb, G4ProductionCuts::GetIndex("e-"));
        if (type == "e+Cut")
            cuts->SetProductionCut(numb, G4ProductionCuts::GetIndex("e+"));
        if (type == "protonCut")
            cuts->SetProductionCut(numb, G4ProductionCuts::GetIndex("proton"));
    }

    if (region != NULL)
        region->SetProductionCuts(cuts);
}
