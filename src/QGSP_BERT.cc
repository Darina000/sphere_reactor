
#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "QGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"

QGSP_BERT::QGSP_BERT(G4int ver)
{
    
    G4double lowlimit = 100.0 * CLHEP::keV;
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit, 100. * CLHEP::TeV);

  defaultCutValue = 1.0*CLHEP::mm;
  SetVerboseLevel(ver);

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver) );
    
    
    
    /*
     // Synchroton Radiation & GN Physics
      RegisterPhysics( new G4EmExtraPhysics(ver) );

      // Decays
      RegisterPhysics( new G4DecayPhysics(ver) );

       // Hadron Elastic scattering
      RegisterPhysics( new G4HadronElasticPhysics(ver) );

      // Hadron Physics
      RegisterPhysics( new G4HadronPhysicsQGSP_BERT(ver));

      // Stopping Physics
      RegisterPhysics( new G4StoppingPhysics(ver) );

      // Ion Physics
      RegisterPhysics( new G4IonPhysics(ver));
      
      // Neutron tracking cut
      RegisterPhysics( new G4NeutronTrackingCut(ver));
     */

}
