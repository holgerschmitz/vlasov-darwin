
#include "diaghelper.h"
#include "potential.h"

void ES_EFieldEnergy::Init (Potential *pField_) {
    t=0;
    
    Energy=0;
    pField = pField_;
    Output.open("ElEnergy.out");
}

void ES_EFieldEnergy::Execute () {

    ScalarField &pGrid = pField->GetEx();
    double *pArr = pGrid.data();
    int S = pGrid.getSize();
    
    Energy=0;
    for (int i=0; i<S; ++i) {
        double E = pArr[i];
        Energy += E*E;
    }
    
    if ( (t%3) == 0)
        Output << t << " " << Energy << endl;
    ++t;

}
