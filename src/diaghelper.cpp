
#include "diaghelper.h"
#include "potential.h"

void ES_EFieldEnergy::Init (Potential *pField_) {
    t=0;
    
    Energy=0;
    pField = pField_;
    Output.open("ElEnergy.out");
}

void ES_EFieldEnergy::Execute () {

    ScalarField &Field = pField->GetEx();
    const int *LBound = Field.getLow();
    const int *UBound = Field.getHigh();
    
    Energy=0;
    for (int i=LBound[0]+1; i<UBound[0]; ++i) 
        for (int j=LBound[1]+1; j<UBound[1]; ++j) {
            double E = Field(i,j);
            Energy += E*E;
        }
    
//    if ( (t%3) == 0)
        Output << t << " " << Energy << endl;
    ++t;

}
