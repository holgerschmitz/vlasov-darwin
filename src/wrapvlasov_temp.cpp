// -*- C++ -*-
// $Id$

#include "vlasovinit.h"
#include "potential.h"
#include "darwin.h"
#include "globals.h"
#include "process.h"

// ----------------------------------------------------------------------------
// VlasovSpecies
// ----------------------------------------------------------------------------
//constructor
//sets densityGoal and dt
template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
VlasovSpecies<ForceField,Advancer,Scheme>::VlasovSpecies(SpeciesData &data)
    : Advancer<ForceField,Scheme>(data)
{
  dt = Parameters::instance().dt();
  densityGoal = data.densityGoal;
}
//-----------------------------------------------------------------------------
//initializer
template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::Init() {
  //reset f_infty, to get interpolationinit called
  this->f_infty = 0;
  
  std::cout << "Initializing: Vlasov Species ..." << std::endl;
  //seems redundant
  PhasePositionI Size;
  //get boundaries
  PositionI Lowx = Process::instance().getBoundary().scalarLow();
  PositionI Highx = Process::instance().getBoundary().scalarHigh();
    
  std::cout << "                             *** RESIZING ***" << std::endl;
  
  //resize
  resize(this->boundary->DistLow(),this->boundary->DistHigh());
  //initialise forcefield
  ForceFieldBase::initialise();
  
  //get potential field and add this
  this->pPot = Parameters::instance().getField();
  this->pPot->AddSpecies(this);

  std::cout << "Resizing fields" << std::endl;

  //resize the kinetic energy matrix
  EKin.resize(Lowx.Data(),Highx.Data());
  
  //get the global grid spacing
  this->dx[0] = Parameters::instance().gridSpace_x();
  this->dx[1] = Parameters::instance().gridSpace_y();
 
  //set the local boxrange accordingly
  for (int i=0; i<2; ++i) this->BoxRange[i]=(Highx[i]-Lowx[i]-1)*this->dx[i];

  std::cout << "Init Base " << this->BoxRange << " " << this->VRange << std::endl;

  std::cout << "Init Field " << this->Charge/this->Mass << " " << this->VRange << std::endl;

  //call init, which calls interpolationinit, which in turn will reset f_infty
  ForceField::Init(this->dx[0]/dt);

  //retrieve the fields from the DiagList, set the fields accordingly
  for (
    ForceFieldBase::DerivedDiagList::iterator it=this->diaglist.begin();
    it !=this->diaglist.end();
    ++it
  )
  {
    (*it)->retrieveField(this);
  }
}
//-----------------------------------------------------------------------------
//wirtes distribution function to stream
template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::write(ostream &O) {
    //get limits
    const int *L = this->Distribution.getLow();
    const int *H = this->Distribution.getHigh();
    //write each matrix element to the stream
    for (int i=L[0]; i<=H[0]; ++i)
      for (int j=L[1]; j<=H[1]; ++j)
        for (int k=L[2]; k<=H[2]; ++k)
          for (int l=L[3]; l<=H[3]; ++l) {
            for (int m=L[4]; m<=H[4]; ++m) 
                O << i << " " << j << " " << 
                    k << " " << l << " " << 
                    m << " " << this->Distribution(i,j,k,l,m) << endl;
                O << endl;
          }

}
//------------------------------------------------------------------------------
//gets called if at calling execute f_infty is zero
template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::InterpolationInitStep(const VlasovDist &Dist) {
  //seems to be redundant, since this gets called if f_infty == 0, this would hold automatically
  this->f_infty = 0;
  //get limits
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
  //set f_infty to the maximum of all distribution values
  for (int i=LBound[0]; i<=UBound[0]; ++i)
    for (int j=LBound[1]; j<=UBound[1]; ++j)
      for (int k=LBound[2]; k<=UBound[2]; ++k) 
        for (int l=LBound[3]; l<=UBound[3]; ++l) 
          for (int m=LBound[4]; m<=UBound[4]; ++m) 
            this->f_infty = max(this->f_infty,this->Distribution(i,j,k,l,m));
  //collect the f_infty value from all processes and set f_infty to the maximum of all results
  this->f_infty = this->boundary->MaxReduce(this->f_infty);
}
//-----------------------------------------------------------------------------
//  Calls the advance method of the Advancer base class.
// then calls Task::Execute to execute sub-Tasks
 
template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::Execute () {
//    if ( (this->tstep%20) == 0 ) {
//        double err = densityError();
//        if (densityGoal!=0) {
//          double err = this->densityError(this->Distribution);
//          this->correctDensityError(err, this->Distribution);
//        } 
//    }
    this->tstep++;
    //if f_infty == 0, then call InitStep, since this is the first call of Init
    if (this->f_infty == 0) InterpolationInitStep(this->Distribution);
//    cerr << "Advance\n";
    Advancer<ForceField,Scheme>::advance(dt);
//    cerr << "Derived Fields\n";
    this->derivedFields.update(*this);
}
//-----------------------------------------------------------------------------


