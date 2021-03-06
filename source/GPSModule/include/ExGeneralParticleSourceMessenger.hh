//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       ExGeneralParticleSourceMessenger.hh
//
// Version:      2.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 2.0, 05/02/2004, Fan Lei, Created.
//    After changes to version 1.1 as in Geant4 v6.0
//     - Mutilple particle source definition
//     - Re-structured commands
//     - old commonds have been retained for backward compatibility, will be
//       removed in the future. 
//
// 2016, Linyan, Modified
// To create UI for correlated events as well background radioactivities
// See documentation - Generator for grammar
//
///////////////////////////////////////////////////////////////////////////////
//
//
// Class Description:
//
// The function of the ExGeneralParticleSourceMessenger is to allow the user to
// enter commands either in interactive command line mode or through macros to
// control the ExGeneralParticleSource. 
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// ExGeneralParticleSourceMessenger(ExGeneralParticleSource *fPtclGun)
//     Constructor:  Sets up commands.
//
// ~ExGeneralParticleSourceMessenger()
//     Destructor:  Deletes commands.
//
// void SetParticleGun(ExSingleParticleSource *fpg) { fParticleGun = fpg; } ;
//     To selecte the particle gun to be defined/modified. 
// void SetNewValue(G4UIcommand *command, G4String newValues)
//     Uses the appropriate methods in the ExGeneralParticleSource to carry out
//     the user commands.
// G4String GetCurrentValue(G4UIcommand *command)
//     Allows the user to retrieve the current values of parameters.
//     Not implemented yet.
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef ExGeneralParticleSourceMessenger_h
#define ExGeneralParticleSourceMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class ExSingleParticleSource;
class ExGeneralParticleSource;

class ExGeneralParticleSourceMessenger: public G4UImessenger
{
  public:
    static ExGeneralParticleSourceMessenger *GetInstance()
    {
        if (ExGeneralParticleSourceMessenger::fInstance == NULL)
            ExGeneralParticleSourceMessenger::fInstance = new ExGeneralParticleSourceMessenger();
        return ExGeneralParticleSourceMessenger::fInstance;
    }
    static ExGeneralParticleSourceMessenger *GetInstance(ExGeneralParticleSource* fExGPS)
    {
        if (ExGeneralParticleSourceMessenger::fInstance == NULL)
            ExGeneralParticleSourceMessenger::fInstance = new ExGeneralParticleSourceMessenger(fExGPS);
        return ExGeneralParticleSourceMessenger::fInstance;
    }

    ExGeneralParticleSourceMessenger(ExGeneralParticleSource* fGPS = NULL);
    ~ExGeneralParticleSourceMessenger();

    void SetParticleGun(ExSingleParticleSource *fpg) { fParticleGun = fpg; } ;
    void SetExGPS(ExGeneralParticleSource *val) { fGPS = val; } ;
    // Select the particle gun to be defined/modified
   
    void SetNewValue(G4UIcommand *command, G4String newValues);
    // Identifies the command which has been invoked by the user, extracts the
    // parameters associated with that command (held in newValues), and uses
    // these values with the appropriate member function of ExGeneralParticleSource.

    G4String GetCurrentValue(G4UIcommand *command);

  private:
    void IonCommand(G4String newValues);
    void IonLvlCommand(G4String newValues);

  private:
    static ExGeneralParticleSourceMessenger* fInstance;
    ExGeneralParticleSource* fGPS;
    ExSingleParticleSource* fParticleGun;
    G4ParticleTable* particleTable;
    G4String histtype;
    
  private: //commands
    G4UIdirectory* gpsDirectory;

    // multiple source control commands
    G4UIdirectory              *sourceDirectory;
    G4UIcmdWithADouble         *addsourceCmd;
    G4UIcmdWithoutParameter    *listsourceCmd;
    G4UIcmdWithoutParameter    *clearsourceCmd;
    G4UIcmdWithoutParameter    *getsourceCmd;
    G4UIcmdWithAnInteger       *setsourceCmd;  
    G4UIcmdWithADouble         *setintensityCmd;
    G4UIcmdWithAnInteger       *deletesourceCmd;
    G4UIcmdWithABool           *multiplevertexCmd;
    G4UIcmdWithABool           *flatsamplingCmd;

    G4UIdirectory              *eventDirectory;
    G4UIcmdWithADouble         *addeventCmd;

    G4UIcmdWithAString         *cgCmd;

    // positional commands
    G4UIdirectory              *positionDirectory;
    G4UIcmdWithAString         *typeCmd1;
    G4UIcmdWithAString         *shapeCmd1;
    G4UIcmdWith3VectorAndUnit  *centreCmd1;
    G4UIcmdWith3Vector         *posrot1Cmd1;
    G4UIcmdWith3Vector         *posrot2Cmd1;
    G4UIcmdWithADoubleAndUnit  *halfxCmd1;
    G4UIcmdWithADoubleAndUnit  *halfyCmd1;
    G4UIcmdWithADoubleAndUnit  *halfzCmd1;
    G4UIcmdWithADoubleAndUnit  *radiusCmd1;
    G4UIcmdWithADoubleAndUnit  *radius0Cmd1;
    G4UIcmdWithADoubleAndUnit  *possigmarCmd1;
    G4UIcmdWithADoubleAndUnit  *possigmaxCmd1;
    G4UIcmdWithADoubleAndUnit  *possigmayCmd1;
    G4UIcmdWithADoubleAndUnit  *paralpCmd1;
    G4UIcmdWithADoubleAndUnit  *partheCmd1;
    G4UIcmdWithADoubleAndUnit  *parphiCmd1;  
    G4UIcmdWithAString         *confineCmd1;
         
  //old ones, will be reomved soon
  G4UIcmdWithAString         *typeCmd;
  G4UIcmdWithAString         *shapeCmd;
  G4UIcmdWith3VectorAndUnit  *centreCmd;
  G4UIcmdWith3Vector         *posrot1Cmd;
  G4UIcmdWith3Vector         *posrot2Cmd;
  G4UIcmdWithADoubleAndUnit  *halfxCmd;
  G4UIcmdWithADoubleAndUnit  *halfyCmd;
  G4UIcmdWithADoubleAndUnit  *halfzCmd;
  G4UIcmdWithADoubleAndUnit  *radiusCmd;
  G4UIcmdWithADoubleAndUnit  *radius0Cmd;
  G4UIcmdWithADoubleAndUnit  *possigmarCmd;
  G4UIcmdWithADoubleAndUnit  *possigmaxCmd;
  G4UIcmdWithADoubleAndUnit  *possigmayCmd;
  G4UIcmdWithADoubleAndUnit  *paralpCmd;
  G4UIcmdWithADoubleAndUnit  *partheCmd;
  G4UIcmdWithADoubleAndUnit  *parphiCmd;  
  G4UIcmdWithAString         *confineCmd; 
        
    // angular commands
    G4UIdirectory* angularDirectory;
    G4UIcmdWithAString         *angtypeCmd1;
    G4UIcmdWith3Vector         *angrot1Cmd1;
    G4UIcmdWith3Vector         *angrot2Cmd1;
    G4UIcmdWithADoubleAndUnit  *minthetaCmd1;
    G4UIcmdWithADoubleAndUnit  *maxthetaCmd1;
    G4UIcmdWithADoubleAndUnit  *minphiCmd1;
    G4UIcmdWithADoubleAndUnit  *maxphiCmd1;
    G4UIcmdWithADoubleAndUnit  *angsigmarCmd1;
    G4UIcmdWithADoubleAndUnit  *angsigmaxCmd1;
    G4UIcmdWithADoubleAndUnit  *angsigmayCmd1;
    G4UIcmdWith3VectorAndUnit  *angfocusCmd;
    G4UIcmdWithABool           *useuserangaxisCmd1;
    G4UIcmdWithABool           *surfnormCmd1;

  // old ones, will be removed soon
  G4UIcmdWithAString         *angtypeCmd;
  G4UIcmdWith3Vector         *angrot1Cmd;
  G4UIcmdWith3Vector         *angrot2Cmd;
  G4UIcmdWithADoubleAndUnit  *minthetaCmd;
  G4UIcmdWithADoubleAndUnit  *maxthetaCmd;
  G4UIcmdWithADoubleAndUnit  *minphiCmd;
  G4UIcmdWithADoubleAndUnit  *maxphiCmd;
  G4UIcmdWithADoubleAndUnit  *angsigmarCmd;
  G4UIcmdWithADoubleAndUnit  *angsigmaxCmd;
  G4UIcmdWithADoubleAndUnit  *angsigmayCmd;
  G4UIcmdWithABool           *useuserangaxisCmd;
  G4UIcmdWithABool           *surfnormCmd;

  // time commands
    G4UIdirectory* timeDirectory;
    G4UIcmdWithAString         *timetypeCmd1;
    G4UIcmdWithADoubleAndUnit  *tminCmd1;
    G4UIcmdWithADoubleAndUnit  *tmaxCmd1;
    G4UIcmdWithADoubleAndUnit  *trangeCmd1;
    G4UIcmdWithADoubleAndUnit  *monotimeCmd1;
    G4UIcmdWithADoubleAndUnit  *tigsigmaCmd1;
    G4UIcmdWithABool           *timespecCmd1;
    G4UIcmdWithADoubleAndUnit         *talphaCmd1;
    G4UIcmdWithADoubleAndUnit         *tzeroCmd1;
    G4UIcmdWithADouble         *tinterceptCmd1;
    G4UIdirectory* timestaDirectory;
	 G4UIcmdWithADouble		 *tssecCmd1;
	 G4UIcmdWithADouble		 *tsnsecCmd1;
    G4UIdirectory* timeendDirectory;
	 G4UIcmdWithADouble		 *tesecCmd1;
	 G4UIcmdWithADouble		 *tensecCmd1;

    // energy commands
    G4UIdirectory* energyDirectory;
    G4UIcmdWithAString         *energytypeCmd1;
    G4UIcmdWithAString         *energyFileCmd;
    G4UIcmdWithAString         *energyHistCmd;
    G4UIcmdWithAString         *energyHistUnitCmd;
    G4UIcmdWithADoubleAndUnit  *eminCmd1;
    G4UIcmdWithADoubleAndUnit  *emaxCmd1;
    G4UIcmdWithADoubleAndUnit  *monoenergyCmd1;
    G4UIcmdWithADoubleAndUnit  *engsigmaCmd1;
    G4UIcmdWithADouble         *alphaCmd1;
    G4UIcmdWithADouble         *tempCmd1;
    G4UIcmdWithADouble         *ezeroCmd1;
    G4UIcmdWithADouble         *gradientCmd1;
    G4UIcmdWithADouble         *interceptCmd1;
    G4UIcmdWithADouble         *arbeintCmd1;
    G4UIcmdWithoutParameter    *calculateCmd1;
    G4UIcmdWithABool           *energyspecCmd1;
    G4UIcmdWithABool           *diffspecCmd1;

    // old ones, will be removed soon
    G4UIcmdWithAString         *energytypeCmd;
    G4UIcmdWithADoubleAndUnit  *eminCmd;
    G4UIcmdWithADoubleAndUnit  *emaxCmd;
    G4UIcmdWithADoubleAndUnit  *monoenergyCmd;
    G4UIcmdWithADoubleAndUnit  *engsigmaCmd;
    G4UIcmdWithADouble         *alphaCmd;
    G4UIcmdWithADouble         *tempCmd;
    G4UIcmdWithADouble         *ezeroCmd;
    G4UIcmdWithADouble         *gradientCmd;
    G4UIcmdWithADouble         *interceptCmd;
    G4UIcmdWithoutParameter    *calculateCmd;
    G4UIcmdWithABool           *energyspecCmd;
    G4UIcmdWithABool           *diffspecCmd;

    // histogram commands
    G4UIdirectory              *histDirectory;
    G4UIcmdWith3Vector         *histpointCmd;
    G4UIcmdWithAString         *histnameCmd;
    G4UIcmdWithAString         *arbintCmd;
    G4UIcmdWithAString         *resethistCmd;
    // old ones, will be removed soon
    G4UIcmdWith3Vector         *histpointCmd1;
    G4UIcmdWithAString         *histfileCmd1;
    G4UIcmdWithAString         *histnameCmd1;
    G4UIcmdWithAString         *arbintCmd1;
    G4UIcmdWithAString         *resethistCmd1;

    G4UIcmdWithAnInteger* verbosityCmd;

    // Commands from G4ParticleGun
    G4UIcommand* ionCmd;
    G4UIcommand* ionLvlCmd;
    G4UIcmdWithAString* radCmd;
    G4UIcmdWithAString* raddetailCmd;
    G4UIcmdWithAString* particleCmd;
    G4UIcmdWithADoubleAndUnit* timeCmd;
    G4UIcmdWith3Vector* polCmd;
    G4UIcmdWithAnInteger* numberCmd;
    G4UIcmdWith3VectorAndUnit* positionCmd;
    G4UIcmdWith3Vector* directionCmd;
    G4UIcmdWithADoubleAndUnit* energyCmd;
    G4UIcmdWithoutParameter* listCmd;

  private: // for ion shooting
    G4bool   fShootIon; 
    G4int    fAtomicNumber;
    G4int    fAtomicMass;
    G4int    fIonCharge;
    G4double fIonExciteEnergy;

    G4bool   fShootIonL;
    G4int    fAtomicNumberL;
    G4int    fAtomicMassL;
    G4int    fIonChargeL;
    G4int    fIonEnergyLevel;
};

#endif

