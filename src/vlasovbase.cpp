#include "vlasovbase.h"
#include "vlasovinit.h"


ForceFieldBase::ForceFieldBase(SpeciesData &data) 
    :  boundary(BoundaryKeeper::getBoundary()),
       tstep(0) { 
    Mass = data.mass;
    Charge = data.charge;
    
    VRange[0] = data.GridRange_vx;
    VRange[1] = data.GridRange_vy;
    VRange[2] = data.GridRange_vz;
    
    init = data.init;
    phasediag = data.phasediag;
}


ForceFieldBase::~ForceFieldBase() {
    std::cout << "Destructing Vlasov Species\n";
    delete boundary;
}


void ForceFieldBase::initialise() {
    std::cout << "Executing initializer " << std::endl;
    
    init->initialise(this);
    
    std::cout << "Boundary Exchange " << std::endl;
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);

    if (NULL != phasediag) phasediag->setField(&Distribution);

    std::cout << "DONE INITIALIZING Vlasov Species" << std::endl;    
}


void ForceFieldBase::resize(PhasePositionI low, PhasePositionI high) {
    std::cout << "RISIZING Vlasov Advancer Base " << low << " " << high << std::endl;
    
    Distribution.resize(low.Data(),high.Data());

    int Pi = 0;
    
    for (int pi = 0; pi<GridSize.length(); ++pi, ++Pi) { 
      GridSize[pi] = high[Pi]-low[Pi]-3;
      GridSizeI[pi] = high[Pi]-low[Pi]-3;     
    }         
    for (int vi = 0; vi < VelSize.length(); ++vi, ++Pi) 
      VelSize[vi] = high[Pi]-low[Pi]; 
                
    VelSizeH = (VelSize-1.0) / 2.0;

    std::cout << "RISIZING Finite Volume Advancer" << std::endl;
    std::cout << "VelSize = " << VelSize << "  VelSizeH = " << VelSizeH << "\n";
}
