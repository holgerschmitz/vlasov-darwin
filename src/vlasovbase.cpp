#include "vlasovbase.h"
#include "vlasovinit.h"


ForceFieldBase::ForceFieldBase(SpeciesData &data) 
    :  boundary(BoundaryKeeper::getBoundary()),
       derivedFields(BoundaryKeeper::getBoundary()),
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
      VelSize[vi] = high[Pi]-low[Pi]-1; 
                
    VelSizeH = (VelSize+1.0) / 2.0;

    std::cout << "RISIZING Finite Volume Advancer" << std::endl;
    std::cout << "VelSize = " << VelSize << "  VelSizeH = " << VelSizeH << "\n";
}


pDistributionDerivedField ForceFieldBase::getDerivedField(std::string name) {
  return derivedFields.getField(name);
}

void ForceFieldBase::addDerivedDiagnostic(VlasovDerivedDiagnostic *ddg) {
  diaglist.push_back(ddg);
}


VlasovDerivedDiagnostic::DerivedDiagList VlasovDerivedDiagnostic::diaglist;
VlasovDerivedDiagnostic *VlasovDerivedDiagnostic::fielddiag=NULL;

VlasovDerivedDiagnostic::VlasovDerivedDiagnostic() {
  diaglist.push_back(this);
}

void VlasovDerivedDiagnostic::retrieveField(ForceFieldBase *vlasov) {
  pDistributionDerivedField dfield = vlasov->getDerivedField(classname);
  if (dfield.pObj()==NULL) 
    setField(NULL);
  else
    setField(&(dfield->getField(fieldname)));
}

PARAMETERMAP* VlasovDerivedDiagnostic::MakeParamMap (PARAMETERMAP* pm) {
  pm = SimpleDiagnostic<ScalarField,std::ofstream>::MakeParamMap(pm);
  (*pm)["class"] = WParameter(new ParameterValue<std::string>(&classname, ""));
  (*pm)["field"] = WParameter(new ParameterValue<std::string>(&fieldname, ""));
  return pm;
}


void CheckDensity(VlasovDist &dist, const char *Msg)
{
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  double Rho = 0;
  
  // The region of the electromagnetic fields is on
  // grid cell smaller in every direction
  // It's important that the density is only calculated on the
  // inner cells, since ScalarFieldReduce simply adds all the densities of
  // the processes.
  
  for (int ix=L[0]+2; ix<=H[0]-2; ++ix) {
    for (int iy=L[1]+2; iy<=H[1]-2; ++iy) {
      for (int j=L[2]; j<=H[2]; ++j)
        for (int k=L[3]; k<=H[3]; ++k)
          for (int l=L[4]; l<=H[4]; ++l) 
              Rho += dist(ix,iy,j,k,l);
    }
  }
  std::cerr << "CheckDensity " << Msg << " " << Rho << std::endl;

}


