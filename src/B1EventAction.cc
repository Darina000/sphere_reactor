
#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{} 

B1EventAction::~B1EventAction()
{}

void B1EventAction::BeginOfEventAction(const G4Event* evt)
{    
  fEdep = 0.;
  G4int evtId = evt->GetEventID();
    //std::cout << "EventID: " << evtId << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);
}
