// -*- C++ -*-
// $Id$

#include "vlasovinit.h"
#include "potential.h"
#include "darwin.h"

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
  dt = Gl_dt;
}

template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::setForceField(typename ForceField::FieldType *pField_) {
  pPot = pField_;
  pPot->AddSpecies(this);
}

template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::Init() {


    cerr << "INITIALIZING Vlasov Species" << endl;

    PhasePositionI Size;
    PositionI Lowx = pPot->GetLBound();
    PositionI Highx = pPot->GetHBound();
    
    cerr << "Resizing" << endl;
    
    resize(boundary->DistLow(),boundary->DistHigh());
    
    cerr << "Resizing fields" << endl;
    
	EKin.resize(Lowx.Data(),Highx.Data());


	dx[0] = GlGridSpace_x;
    dx[1] = GlGridSpace_y;
    
    for (int i=0; i<2; ++i) BoxRange[i]=(Highx[i]-Lowx[i]-3)*dx[i];
    
    cerr << "Init Base " << BoxRange << " " << VRange << endl;
       
       
    cerr << "Init Field " << Charge/Mass << " " << VRange << endl;
    ForceField::Init(dx[0]/dt);
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
void VlasovSpecies<ForceField,Advancer,Scheme>::writeYVySlice(int x,int vx,int vz, string fname) {
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();
    if (   (x<L[0]+2) || (x>H[0]-2) 
        || (vx<L[2]) || (vx>H[2])
        || (vz<L[4]) || (vz>H[4]))  return;
    
    ofstream Slice(fname.c_str());
    
    for (int j=L[1]+2; j<=H[1]-2; ++j)
        for (int k=L[3]; k<=H[3]; ++k)
            Slice << j << " " << k << " " 
                << Distribution(x,j,vx,k,vz) << endl;
    
    Slice.close();
}

template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::writeVxVySlice(int t, int x, int y, int vz, ostream &Slice) 
{
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();
    if (   (x<L[0]+2) || (x>H[0]-2) 
        || (y<L[1]+2) || (y>H[1]-2)
        || (vz<L[4]) || (vz>H[4]))  return;
    
    
    for (int j=L[2]; j<=H[2]; ++j) 
    {
        for (int k=L[3]; k<=H[3]; ++k)
        {
            Slice << t << " " << j << " " << k << " " 
                << Distribution(x,y,j,k,vz) << endl;
        }
        Slice << endl;
    }
          
    
}



template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
ScalarField &VlasovSpecies<ForceField,Advancer,Scheme>::KineticEnergy () {
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();

    VelocityI vi;
    VelocityD V;
    double d;

    EKin.clear();

    for (int i=L[0]+2; i<=H[0]-2; ++i) 
      for (int j=L[1]+2; j<=H[1]-2; ++j)

        for (vi[0]=L[2]; vi[0]<=H[2]; ++vi[0])
          for (vi[1]=L[3]; vi[1]<=H[3]; ++vi[1])
            for (vi[2]=L[4]; vi[2]<=H[4]; ++vi[2]) {
                V = velocity(vi);
                d = Distribution(i,j,vi[0],vi[1],vi[2]);
                EKin(i,j) += (V[0]*V[0]+V[1]*V[1]+V[2]*V[2])*d;
    }
    
    boundary->ScalarFieldReduce(EKin);
    return EKin;
}

template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
double VlasovSpecies<ForceField,Advancer,Scheme>::TotalEnergy () {
    ScalarField &kin = KineticEnergy();
    ScalarField &fld = FieldEnergy();
    const int *L = kin.getLow();
    const int *H = kin.getHigh();
    
    double ETot = 0;
    double dxdy = dx[0]*dx[1];
    
    
    for (int i=L[0]+1; i<=H[0]-1; ++i) 
        for (int j=L[1]+1; j<=H[1]-1; ++j) {
            ETot += kin(i,j)+ fld(i,j);
    }
    return ETot;
}


template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::InterpolationInitStep(const VlasovDist &Dist) {
    f_infty = 1;
}


template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
double VlasovSpecies<ForceField,Advancer,Scheme>::densityError() {
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    double avg = 0;

    for (int i=LBound[0]+2; i<=UBound[0]-2; ++i)
      for (int j=LBound[1]+2; j<=UBound[1]-2; ++j)
        for (int k=LBound[2]; k<=UBound[2]; ++k) 
          for (int l=LBound[3]; l<=UBound[3]; ++l) 
            for (int m=LBound[4]; m<=UBound[4]; ++m) 
              avg += Distribution(i,j,k,l,m);
    
    avg = avg/double((UBound[0]-LBound[0]-3)*(UBound[1]-LBound[1]-3));
    
    cerr << "Partial Error in density: " << avg - 1 << endl;

    avg =  boundary->AvgReduce(avg);
    cerr << "Total Error in density: " << avg - 1 << endl;
    return avg;
}

template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
void VlasovSpecies<ForceField,Advancer,Scheme>::correctDensityError(double err) {
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    for (int i=LBound[0]; i<=UBound[0]; ++i)
      for (int j=LBound[1]; j<=UBound[1]; ++j)
        for (int k=LBound[2]; k<=UBound[2]; ++k) 
          for (int l=LBound[3]; l<=UBound[3]; ++l) 
            for (int m=LBound[4]; m<=UBound[4]; ++m) 
              Distribution(i,j,k,l,m) /= err;
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
    if ( (tstep%20) == 0) {
        double err = densityError();
        correctDensityError(err);
    }
    tstep++;
    InterpolationInitStep(Distribution);
    advance(dt);
    derivedFields.update(*this);
}

template<
  class ForceField, 
  template<class, template<class> class> class Advancer,
  template<class> class Scheme
>
PARAMETERMAP* VlasovSpecies<ForceField,Advancer,Scheme>::MakeParamMap (PARAMETERMAP* pm)
{
  pm = Task::MakeParamMap(pm);
  return Advancer<ForceField,Scheme>::MakeParamMap(pm);
}
