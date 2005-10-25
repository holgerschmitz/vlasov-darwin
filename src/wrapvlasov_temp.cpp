// -*- C++ -*-
// $Id$

#include "vlasovinit.h"
#include "potential.h"
#include "darwin.h"
#include "globals.h"
#include "process.h"

// -------------------------------------------------------------------
// VlasovSpecies
// -------------------------------------------------------------------


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


template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::Init() {

  this->f_infty = 0;
  
  std::cout << "INITIALIZING Vlasov Species" << std::endl;

  PhasePositionI Size;
  PositionI Lowx = Process::instance().getBoundary().scalarLow();
  PositionI Highx = Process::instance().getBoundary().scalarHigh();
    
  std::cout << "Resizing" << std::endl;
    
  resize(this->boundary->DistLow(),this->boundary->DistHigh());
  ForceFieldBase::initialise();
    
  this->pPot = Parameters::instance().getField();
  this->pPot->AddSpecies(this);

  std::cout << "Resizing fields" << std::endl;
    
	EKin.resize(Lowx.Data(),Highx.Data());


	this->dx[0] = Parameters::instance().gridSpace_x();
  this->dx[1] = Parameters::instance().gridSpace_y();
    
  for (int i=0; i<2; ++i) this->BoxRange[i]=(Highx[i]-Lowx[i]-1)*this->dx[i];
    
  std::cout << "Init Base " << this->BoxRange << " " << this->VRange << std::endl;
       
       
  std::cout << "Init Field " << this->Charge/this->Mass << " " << this->VRange << std::endl;
  ForceField::Init(this->dx[0]/dt);
  
  for (
    ForceFieldBase::DerivedDiagList::iterator it=this->diaglist.begin();
    it !=this->diaglist.end();
    ++it
  )
  {
    (*it)->retrieveField(this);
  }
}



template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::write(ostream &O) {
    const int *L = this->Distribution.getLow();
    const int *H = this->Distribution.getHigh();

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



template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::InterpolationInitStep(const VlasovDist &Dist) {
  this->f_infty = 0;
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
  
  for (int i=LBound[0]; i<=UBound[0]; ++i)
    for (int j=LBound[1]; j<=UBound[1]; ++j)
      for (int k=LBound[2]; k<=UBound[2]; ++k) 
        for (int l=LBound[3]; l<=UBound[3]; ++l) 
          for (int m=LBound[4]; m<=UBound[4]; ++m) 
            this->f_infty = max(this->f_infty,this->Distribution(i,j,k,l,m));
            
  this->f_infty = this->boundary->MaxReduce(this->f_infty);
}





/** @brief Calls the advance method of the Advancer base class.
 * then calls Task::Execute to execute sub-Tasks
 */
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
    if (this->f_infty == 0) InterpolationInitStep(this->Distribution);
//    cerr << "Advance\n";
    Advancer<ForceField,Scheme>::advance(dt);
//    cerr << "Derived Fields\n";
    this->derivedFields.update(*this);
}



