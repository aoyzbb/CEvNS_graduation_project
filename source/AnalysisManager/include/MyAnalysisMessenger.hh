//*********************************************
//  This is Geant4 Template
//                                  author:Qian
//

#ifndef MyAnalysisMessenger_h
#define MyAnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class MyRunAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class MyAnalysisMessenger : public G4UImessenger
{
public:
    MyAnalysisMessenger();
    ~MyAnalysisMessenger();

    virtual void SetNewValue(G4UIcommand *, G4String);

private:
    G4UIcmdWithAString *fFileNameCmd1;
    G4UIcmdWithAString *fFileNameCmd2;
    G4UIcmdWithoutParameter *fSelectCmd1;
    G4UIcmdWithoutParameter *fSelectCmd2;
    G4UIcmdWithoutParameter *fSelectCmd3;
    G4UIcmdWithoutParameter *fSelectCmd4;
};

#endif
