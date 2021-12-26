#include "MyRunManager.hh"

MyRunManager::MyRunManager(): G4RunManager(){
    fPGGeneratorList = PGGeneratorList::GetInstance();
}
MyRunManager::~MyRunManager(){}

void MyRunManager::GenEvents(G4int GenEvents, G4int GenValidEvents)
{
    if(GenEvents > 0)
        this->BeamOn(GenEvents);
    else
    {
        this->GenValidEvents(GenValidEvents);
    }
}

void MyRunManager::GenValidEvents(G4int GenValidEvents)
{
    //runManager->InitializeEventLoop(n_event,macroFile,n_select);
    G4int ValidEvents = 0;
    this->ConstructScoringWorlds();
    this->RunInitialization();
    for(G4int i_event = 0; ValidEvents < GenValidEvents; ++i_event)
    {
        this->ProcessOneEvent(i_event);
        this->TerminateOneEvent();
        if(fPGGeneratorList->IfValidEvent()) 
            ++ValidEvents;
    }
    this->TerminateEventLoop();
    this->RunTermination();
}