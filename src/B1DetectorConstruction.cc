
#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include <G4UniformMagField.hh>
#include <G4ClassicalRK4.hh>
#include <G4TransportationManager.hh>
#include <G4EqMagElectricField.hh>
#include "G4EquationOfMotion.hh"

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"
#include "G4ElectroMagneticField.hh"
#include "G4MonopoleEq.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4EqEMFieldWithSpin.hh"
#include "G4PropagatorInField.hh"

F05Field::F05Field() : G4ElectroMagneticField()
{
}

F05Field::~F05Field()
{
}

void F05Field::GetFieldValue( const G4double Point[4], G4double* Bfield ) const
{
  // Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time

  const G4double Bz = 0.000*tesla;
  const G4double Er = 520.000*kilovolt/m;

  G4double Ex, Ey, Ez;
    
    double x, y, z;
    x = Point[0]/m;
    y = Point[1]/m;
    z = Point[2]/m;
    
    
  G4double R0 = std::sqrt(std::pow(x,2) + std::pow(y,2)+std::pow(z,2));

    double cos_theta, cos_phi, sin_theta, sin_phi;
    
    if (R0>0 ){
        //coefficient for recalculating the field value for each point in the case of a uniform field
        double k = 1/R0;
        
        Ex = 1*x*k* Er;
        Ey = 1*y*k* Er;
        Ez = 1*z*k* Er;

    }else{
       Ex=0;
       Ey=0;
        Ez=0;
    }
  //magnetic field
  Bfield[0]=0;
  Bfield[1]=0;
  Bfield[2]=0;

    //Electric field
  Bfield[3]=Ex;
  Bfield[4]=Ey;
  Bfield[5]=Ez;

  return;
}


B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

B1DetectorConstruction::~B1DetectorConstruction()
{ }

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
    
    G4String name, symbol;
    G4double density, fractionmass;
    G4int nel, ncomponents;
    
    G4NistManager* nist = G4NistManager::Instance();
    G4int nvar = 8;
    fMinStep = 0.001*mm ;
    G4bool isotopes = true;
    
    G4bool checkOverlaps = true;
    auto nistManager = G4NistManager::Instance();
      
    density = 0.4113*kg/m3;
    auto AirDensity = nistManager ->BuildMaterialWithNewDensity("AirDensity", "G4_AIR", density);

  // World
  //
  G4double world_sizeXY = 250*m;
  G4double world_sizeZ  = 250*m;

    G4Box* solidWorld =
    new G4Box("World",
              1.0*world_sizeXY, 1.0*world_sizeXY, 1.0*world_sizeZ);
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,
                        AirDensity,
                        "World");
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, 0),
                      logicWorld,
                      "World",
                      0,
                      false,
                      0,
                      checkOverlaps);

    
    
//cap
    double a = 250.0;// world_sizeZ*;
    G4Sphere* solidShape2_1 =
      new G4Sphere("Shape2_1", (a-1)*m, a*m, 1, 1, 1, 1);
                  
    G4LogicalVolume* logicShape2_1 =
      new G4LogicalVolume(solidShape2_1,
                          AirDensity,          //its material
                          "Shape2_1");
                 
    new G4PVPlacement(0,
                      G4ThreeVector(0*cm, 0, 0*cm),
                      logicShape2_1,
                      "Shape2_1",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
  
   

    //uniform field in logicWorld
    
    /*
    G4ElectricField* magField_5 = new G4UniformElectricField( G4ThreeVector(0.0,0.0,
                                                                              -150*kilovolt/m));
        G4EqMagElectricField* fEquation1_5 = new G4EqMagElectricField(magField_5);
        G4MagIntegratorStepper* fStepper1_5 = new G4ClassicalRK4( fEquation1_5, nvar );
        G4FieldManager* localFieldMgr1_5
                = new G4FieldManager(magField_5);
    logicWorld->SetFieldManager(localFieldMgr1_5, true);
        G4MagInt_Driver* fIntgrDriver1_5 = new G4MagInt_Driver(fMinStep,
                                         fStepper1_5,
                                         fStepper1_5->GetNumberOfVariables() );
        G4ChordFinder* fChordFinder1_5 = new G4ChordFinder(fIntgrDriver1_5);
           localFieldMgr1_5->SetChordFinder( fChordFinder1_5 );
    */
    
  fScoringVolume = logicWorld;

  return physWorld;
}

G4ThreadLocal F05Field* B1DetectorConstruction::fField = 0;

void B1DetectorConstruction::ConstructSDandField(){
    if (!fField) {
       fField = new F05Field();
       G4EqEMFieldWithSpin* equation = new G4EqEMFieldWithSpin(fField);

       G4FieldManager* fieldManager
        = G4TransportationManager::GetTransportationManager()->GetFieldManager();

       fieldManager->SetDetectorField(fField);

       G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,12);

       G4double minStep           = 0.01*mm;

       G4ChordFinder* chordFinder =
                      new G4ChordFinder((G4MagneticField*)fField,minStep,stepper);

       // Set accuracy parameters
       G4double deltaChord        = 3.0*mm;
       chordFinder->SetDeltaChord( deltaChord );

       G4double deltaOneStep      = 0.001*mm;
       fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

       G4double deltaIntersection = 0.01*mm;
       fieldManager->SetDeltaIntersection(deltaIntersection);

       G4TransportationManager* transportManager =
                             G4TransportationManager::GetTransportationManager();

       G4PropagatorInField* fieldPropagator =
                                        transportManager->GetPropagatorInField();

       G4double epsMin            = 2.5e-5*mm;
       G4double epsMax            = 0.09*mm;

       fieldPropagator->SetMinimumEpsilonStep(epsMin);
       fieldPropagator->SetMaximumEpsilonStep(epsMax);

       fieldManager->SetChordFinder(chordFinder);
        
        
    }
    
}
