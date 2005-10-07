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

  f_infty = 0;
  
  std::cout << "INITIALIZING Vlasov Species" << std::endl;

  PhasePositionI Size;
  PositionI Lowx = Process::instance().getBoundary().scalarLow();
  PositionI Highx = Process::instance().getBoundary().scalarHigh();
    
  std::cout << "Resizing" << std::endl;
    
  resize(boundary->DistLow(),boundary->DistHigh());
  ForceFieldBase::initialise();
    
  pPot = Parameters::instance().getField();
  pPot->AddSpecies(this);

  std::cout << "Resizing fields" << std::endl;
    
	EKin.resize(Lowx.Data(),Highx.Data());


	dx[0] = Parameters::instance().gridSpace_x();
  dx[1] = Parameters::instance().gridSpace_y();
    
  for (int i=0; i<2; ++i) BoxRange[i]=(Highx[i]-Lowx[i]-1)*dx[i];
    
  std::cout << "Init Base " << BoxRange << " " << VRange << std::endl;
       
       
  std::cout << "Init Field " << Charge/Mass << " " << VRange << std::endl;
  ForceField::Init(dx[0]/dt);
  
  for (
    ForceFieldBase::DerivedDiagList::iterator it=diaglist.begin();
    it !=diaglist.end();
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
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();

    for (int i=L[0]; i<=H[0]; ++i)
      for (int j=L[1]; j<=H[1]; ++j)
        for (int k=L[2]; k<=H[2]; ++k)
          for (int l=L[3]; l<=H[3]; ++l) {
            for (int m=L[4]; m<=H[4]; ++m) 
                O << i << " " << j << " " << 
                    k << " " << l << " " << 
                    m << " " << Distribution(i,j,k,l,m) << endl;
            O << endl;
          }
            
}



template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::InterpolationInitStep(const VlasovDist &Dist) {
  f_infty = 0;
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();
  
  for (int i=LBound[0]; i<=UBound[0]; ++i)
    for (int j=LBound[1]; j<=UBound[1]; ++j)
      for (int k=LBound[2]; k<=UBound[2]; ++k) 
        for (int l=LBound[3]; l<=UBound[3]; ++l) 
          for (int m=LBound[4]; m<=UBound[4]; ++m) 
            f_infty = max(f_infty,Distribution(i,j,k,l,m));
            
  f_infty = boundary->MaxReduce(f_infty);
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
    if ( (tstep%20) == 0 ) {
//        double err = densityError();
        if (densityGoal!=0) {
          double err = densityError(Distribution);
          correctDensityError(err, Distribution);
        } 
    }
    tstep++;
    if (f_infty == 0) InterpolationInitStep(Distribution);
//    cerr << "Advance\n";
    Advancer<ForceField,Scheme>::advance(dt);
//    cerr << "Derived Fields\n";
    derivedFields.update(*this);
}



