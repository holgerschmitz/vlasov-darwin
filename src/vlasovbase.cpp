#include "vlasovbase.h"
#include "vlasovinit.h"


ForceFieldBase::ForceFieldBase(SpeciesData &data) 
    :  boundary(data.bound),
       tstep(0) { 
    Mass = data.mass;
    Charge = data.charge;
    
    VRange[0] = data.GridRange_vx;
    VRange[1] = data.GridRange_vy;
    VRange[2] = data.GridRange_vz;
}

PARAMETERMAP* ForceFieldBase::MakeParamMap (PARAMETERMAP* pm) {
  (*pm)["mass"] = WParameter(new ParameterValue<double>(&Mass, 1));
  (*pm)["charge"] = WParameter(new ParameterValue<double>(&Charge, 1));
  (*pm)["single-periodic-boundary"] = WParameter(new ParameterTask<SinglePeriodicBoundary, Boundary>(&boundary));
#ifndef SINGLE_PROCESSOR
  (*pm)["mpi-periodic-split-x-boundary"] = WParameter(new ParameterTask<MPIPeriodicSplitXBoundary, Boundary>(&boundary));
  (*pm)["mpi-periodic-split-xy-boundary"] = WParameter(new ParameterTask<MPIPeriodicSplitXYBoundary, Boundary>(&boundary));
#endif
}

ForceFieldBase::~ForceFieldBase() {
    cerr << "Destructing Vlasov Species\n";
    delete boundary;
}


void ForceFieldBase::initialise(VlasovInitialiser *init) {
    cerr << "Executing initializer " << endl;
    
    init->initialise(Distribution, VRange);
    
    cerr << "Boundary Exchange " << endl;
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);

    cerr << "DONE INITIALIZING Vlasov Species" << endl;
}


void ForceFieldBase::resize(PhasePositionI &low, PhasePositionI &high) {
    cerr << "RISIZING Vlasov Advancer Base " << low << " " << high << endl;
    
    Distribution.resize(low.Data(),high.Data());

    int Pi = 0;
    
    for (int pi = 0; pi<GridSize.length(); ++pi, ++Pi) { 
      GridSize[pi] = high[Pi]-low[Pi]-3;
      GridSizeI[pi] = high[Pi]-low[Pi]-3;     
    }         
    for (int vi = 0; vi < VelSize.length(); ++vi, ++Pi) 
      VelSize[vi] = high[Pi]-low[Pi]; 
                
    VelSizeH = (VelSize-1.0) / 2.0;

    cerr << "RISIZING Finite Volume Advancer" << endl;
    cerr << "VelSize = " << VelSize << "  VelSizeH = " << VelSizeH << "\n";
}
