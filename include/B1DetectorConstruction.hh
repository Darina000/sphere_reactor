
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include <G4MagIntegratorStepper.hh>
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4MagIntegratorDriver.hh>
#include<G4VSolid.hh>
#include<G4LogicalVolume.hh>
#include<G4VPhysicalVolume.hh>
#include<G4Material.hh>

class G4VPhysicalVolume;
class G4LogicalVolume;

class F05Field : public G4ElectroMagneticField
{

public:

  F05Field();
  virtual ~F05Field();

  // DoesFieldChangeEnergy() returns true.
  virtual G4bool DoesFieldChangeEnergy() const { return true; };

  // GetFieldValue() returns the field value at a given point[].
  // field is really field[6]: Bx,By,Bz,Ex,Ey,Ez.
  // point[] is in global coordinates: x,y,z,t.
  virtual void GetFieldValue(const G4double Point[4], G4double* Bfield) const;
};


class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    
    
    G4ElectricField*        fEMfield;
    G4EqMagElectricField*   fEquation;
    G4MagIntegratorStepper* fStepper;
    G4FieldManager*         fFieldMgr;
    G4double                fMinStep ;
    G4ChordFinder*          fChordFinder ;
    G4MagInt_Driver*        fIntgrDriver;

  protected:
    G4LogicalVolume*  fScoringVolume;
private:
    static G4ThreadLocal F05Field* fField;
};
