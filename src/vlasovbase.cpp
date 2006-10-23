#include "vlasovbase.h"
#include "vlasovinit.h"
//-----------------------------------------------------------------------------
//ForceFieldBase
//-----------------------------------------------------------------------------

// perform initialisation of the derived fields from the boundary
// then assign the data members of ForceFeldBase with the values from species data
// set timestep to zero
// todo wouldnt an init list work for the data members?
ForceFieldBase::ForceFieldBase(SpeciesData &data) 
    :  boundary(BoundaryKeeper::getBoundary()),
       derivedFields(BoundaryKeeper::getBoundary()),
       tstep(0) { 
  //get the values from species data
  Mass = data.mass;
  Charge = data.charge;
  
  VRange[0] = data.GridRange_vx;
  VRange[1] = data.GridRange_vy;
  VRange[2] = data.GridRange_vz;
  
  init = data.init;
  phasediag = data.phasediag;
}
//-----------------------------------------------------------------------------
//delete the boundary object from the heap
ForceFieldBase::~ForceFieldBase() {
    std::cout << "Destructing Vlasov Species\n";
    delete boundary;
}

//-----------------------------------------------------------------------------
void ForceFieldBase::initialise() {
    std::cout << "Executing initializer " << std::endl;
    
    //call the pointed to initialiser
    init->initialise(this);
    
    //exchange boundary with neighbouring processes
    std::cout << "Boundary Exchange " << std::endl;
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);

    if (NULL != phasediag) phasediag->setField(&Distribution);

    std::cout << "Initializing Vlasov Species: done" << std::endl;    
}
//-----------------------------------------------------------------------------

void ForceFieldBase::resize(PhasePositionI low, PhasePositionI high) {
	std::cout << "Resizing: Vlasov Advancer Base " << low << " " << high << std::endl;
    //resizes the distribution matrix
    Distribution.resize(low.Data(),high.Data());

    int Pi = 0;
    
    for (int pi = 0; pi<GridSize.length(); ++pi, ++Pi) { 
      GridSize[pi] = high[Pi]-low[Pi]-3;
      GridSizeI[pi] = high[Pi]-low[Pi]-3;     
    }         
    for (int vi = 0; vi < VelSize.length(); ++vi, ++Pi)
    {
      VelSize[vi] = high[Pi]-low[Pi]-1;
      VelGridCenter[vi]  = 0.5*(high[Pi]+low[Pi]);
      VelSizeH[vi] = 0.5*(high[Pi]-low[Pi]);
    }
                

    std::cout << "Resizing: Finite Volume Advancer" << std::endl;
    std::cout << "VelSize = " << VelSize << "  VelSizeH = " << VelSizeH << "\n";
}
//-----------------------------------------------------------------------------

double ForceFieldBase::densityError(VlasovDist &distrib) {

    //get the index limits
    const int *UBound = distrib.getHigh();
    const int *LBound = distrib.getLow();

    double avg = 0;

    //add all values of the distribution matrix on the INNER cells
    for (int i=LBound[0]+2; i<=UBound[0]-2; ++i)
      for (int j=LBound[1]+2; j<=UBound[1]-2; ++j)
        for (int k=LBound[2]; k<=UBound[2]; ++k) 
          for (int l=LBound[3]; l<=UBound[3]; ++l) 
            for (int m=LBound[4]; m<=UBound[4]; ++m)
            {
              avg += distrib(i,j,k,l,m);
            }
    //divide by the extent in cells  to get the average
    avg = avg/double((UBound[0]-LBound[0]-3)*(UBound[1]-LBound[1]-3));
    //call the reduce method for this value
    //
    avg =  boundary->AvgReduce(avg);
    //if this is the master process write the total density
    if (boundary->master())
      std::cout << "Total density: " << avg  << std::endl;
    
    return avg;
}
//-----------------------------------------------------------------------------

void ForceFieldBase::correctDensityError(double err, VlasovDist &distrib) {
    if (boundary->master())
      std::cout << "Correcting density: " << err  << std::endl;
    
    //divide all values by the error
    const int *UBound = distrib.getHigh();
    const int *LBound = distrib.getLow();
    for (int i=LBound[0]; i<=UBound[0]; ++i)
      for (int j=LBound[1]; j<=UBound[1]; ++j)
        for (int k=LBound[2]; k<=UBound[2]; ++k) 
          for (int l=LBound[3]; l<=UBound[3]; ++l) 
            for (int m=LBound[4]; m<=UBound[4]; ++m) 
              distrib(i,j,k,l,m) /= err;
}
//-----------------------------------------------------------------------------

//return the underliing field by its name
pDistributionDerivedField ForceFieldBase::getDerivedField(std::string name) {
  return derivedFields.getField(name);
}
//-----------------------------------------------------------------------------

//register the diagnostic object by adding it to the diagnostics list
void ForceFieldBase::addDerivedDiagnostic(VlasovDerivedDiagnostic *ddg) {
  diaglist.push_back(ddg);
}
//ForceFieldBase
//-----------------------------------------------------------------------------
//VlasovDerivedDiagnostic

VlasovDerivedDiagnostic::DerivedDiagList VlasovDerivedDiagnostic::diaglist;
VlasovDerivedDiagnostic *VlasovDerivedDiagnostic::fielddiag=NULL;

//constructor, register *this with the diaglist
VlasovDerivedDiagnostic::VlasovDerivedDiagnostic() {
  diaglist.push_back(this);
}

// retrieve field by field- and classname from the ForceFieldBase object
// setField is a method of the SimpleDiagnostic base class
// pDistributionDerivedField is a typedef for a wrapped pointer to DistributionDerivedField
// DistributionDerivedField is the base class for the derived field hierachy
void VlasovDerivedDiagnostic::retrieveField(ForceFieldBase *vlasov) {
  pDistributionDerivedField dfield = vlasov->getDerivedField(classname);
  if (dfield.pObj()==NULL) 
    setField(NULL);
  else
    setField(&(dfield->getField(fieldname)));
}

//creates the parameter map for the Rebuild method of the Rebuildable base class
PARAMETERMAP* VlasovDerivedDiagnostic::MakeParamMap (PARAMETERMAP* pm) {
  pm = SimpleDiagnostic<ScalarField,std::ofstream>::MakeParamMap(pm);
  (*pm)["class"] = WParameter(new ParameterValue<std::string>(&classname, ""));
  (*pm)["field"] = WParameter(new ParameterValue<std::string>(&fieldname, ""));
  return pm;
}
//VlasovDerivedDiagnostic
//-----------------------------------------------------------------------------
//CheckDensity

void CheckDensity(VlasovDist &dist, const char *Msg)
{
  // summation limits, retrieved from the Matrix object dist
  // since in multi-processor runs this matrix is split
  // these values need to be fetched from the local submatrix
  // else index troubles would be the consequence
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  double Rho = 0;
  
  // The region of the electromagnetic fields is on
  // grid cell smaller in every direction
  // It's important that the density is only calculated on the
  // inner cells, since ScalarFieldReduce simply adds all the densities of
  // the processes.
  // This simply adds up the density over all local cells of the 5dim grid
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
//-----------------------------------------------------------------------------


