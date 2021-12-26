
#ifndef MyRunManager_h
#define MyRunManager_h 1

#include "G4RunManager.hh"
#include "PGGeneratorList.hh"

class MyRunManager: public G4RunManager
{
public:
    MyRunManager();
    ~MyRunManager();

    void GenEvents(G4int, G4int);

private:
    void GenValidEvents(G4int);
    PGGeneratorList* fPGGeneratorList;
};

#endif