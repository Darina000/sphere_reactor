
#include "G4VUserActionInitialization.hh"

class B1ActionInitialization : public G4VUserActionInitialization
{
  public:
    B1ActionInitialization();
    virtual ~B1ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

    
