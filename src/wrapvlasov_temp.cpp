// -*- C++ -*-
// $Id$

#include "vlasovinit.h"
#include "potential.h"
#include "darwin.h"

// -------------------------------------------------------------------
// VlasovSpecies
// -------------------------------------------------------------------

template<class ForceField>
VlasovSpecies<ForceField>::VlasovSpecies (SpeciesData &data) 
    :  ForceField(),
       RKState(-2),
       boundary(data.bound),
       t(0) { 
    dt = Gl_dt;
    Mass = data.mass;
    Charge = data.charge;

    
    VRange[0] = data.GridRange_vx;
    VRange[1] = data.GridRange_vy;
    VRange[2] = data.GridRange_vz;
}

template<class ForceField>
VlasovSpecies<ForceField>::~VlasovSpecies () {
    cerr << "Destructing Vlasov Species\n";
    delete boundary;
}


template<class ForceField>
void VlasovSpecies<ForceField>::setForceField(typename ForceField::FieldType *pField_) {
  pPot = pField_;
  pPot->AddSpecies(this);
}

template<class ForceField>
void VlasovSpecies<ForceField>::initialise(VlasovInitialiser *init) {
    cerr << "Executing initializer " << endl;
    
    init->initialise(Distribution, VRange);
    
    cerr << "Boundary Exchange " << endl;
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);

    cerr << "DONE INITIALIZING Vlasov Species" << endl;
}

template<class ForceField>
void VlasovSpecies<ForceField>::Init() {


    cerr << "INITIALIZING Vlasov Species" << endl;

    PhasePositionI Size;
    PositionI Lowx = pPot->GetLBound();
    PositionI Highx = pPot->GetHBound();
    
    cerr << "Resizing" << endl;
    
    resize(boundary->DistLow(),boundary->DistHigh());
    
    cerr << "Resizing fields" << endl;
    
	gRho.resize(Lowx.Data(),Highx.Data());
	EKin.resize(Lowx.Data(),Highx.Data());

	Jx.resize(Lowx.Data(),Highx.Data());
	Jy.resize(Lowx.Data(),Highx.Data());
	Jz.resize(Lowx.Data(),Highx.Data());

	Vxx.resize(Lowx.Data(),Highx.Data());
	Vxy.resize(Lowx.Data(),Highx.Data());
	Vxz.resize(Lowx.Data(),Highx.Data());
	Vyy.resize(Lowx.Data(),Highx.Data());
	Vyz.resize(Lowx.Data(),Highx.Data());
	Vzz.resize(Lowx.Data(),Highx.Data());

	dx[0] = GlGridSpace_x;
    dx[1] = GlGridSpace_y;
    
    for (int i=0; i<2; ++i) BoxRange[i]=(Highx[i]-Lowx[i]-3)*dx[i];
    
    cerr << "Init Base " << BoxRange << " " << VRange << endl;
       
       
    cerr << "Init Field " << Charge/Mass << " " << VRange << endl;
    ForceField::Init(dx[0]/dt);
}

template<class ForceField>
void VlasovSpecies<ForceField>::resize(PhasePositionI &low, PhasePositionI &high) {
    cerr << "RISIZING Vlasov Advancer Base " << low << " " << high << endl;
    
    Distribution.resize(low.Data(),high.Data());

    int Pi = 0;
    
    for (   int pi = 0; pi<GridSize.length(); ++pi, ++Pi) { 
        GridSize[pi] = high[Pi]-low[Pi]-3;
        GridSizeI[pi] = high[Pi]-low[Pi]-3;     
    }         
    for (   int vi = 0; vi < VelSize.length(); ++vi, ++Pi) 
                VelSize[vi] = high[Pi]-low[Pi]; 
                
    VelSizeH = (VelSize-1.0) / 2.0;

    cerr << "RISIZING Finite Volume Advancer" << endl;
    cerr << "VelSize = " << VelSize << "  VelSizeH = " << VelSizeH << "\n";
//    TempDist.resize(low.Data(),high.Data());
//    T1.resize(low.Data(),high.Data());
//    T2.resize(low.Data(),high.Data());
}


template<class ForceField>
void VlasovSpecies<ForceField>::write(ostream &O) {
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

template<class ForceField>
void VlasovSpecies<ForceField>::writeYVySlice(int x,int vx,int vz, string fname) {
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

template<class ForceField>
void VlasovSpecies<ForceField>::writeVxVySlice(int t, int x, int y, int vz, ostream &Slice) 
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


/** @brief Create density with the density() method supplied by
 *  VlasovAdvancerBase and copy the result into a grid, as requested
 *  by the Opar architecture.
 */
template<class ForceField>
void VlasovSpecies<ForceField>::MakeRho () {

//    cerr << "VlasovSpecies::MakeRho" << endl;

    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();
    
    gRho.clear();
    
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
                gRho(ix,iy) += Distribution(ix,iy,j,k,l);
//        cerr << "den " << ix << " " << iy << " " << gRho(ix,iy) << endl;
        }
    }
    
    boundary->ScalarFieldReduce(gRho);
}

template<class ForceField>
void VlasovSpecies<ForceField>::MakeJs () {
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();
    
    VelocityI vi;
    VelocityD V;
    double d;
    
    Jx.clear();
    Jy.clear();
    Jz.clear();

    Vxx.clear();
    Vxy.clear();
    Vxz.clear();
    Vyy.clear();
    Vyz.clear();
    Vzz.clear();

    for (int i=L[0]+2; i<=H[0]-2; ++i) {
      for (int j=L[1]+2; j<=H[1]-2; ++j) {

        for (vi[0]=L[2]; vi[0]<=H[2]; ++vi[0])
          for (vi[1]=L[3]; vi[1]<=H[3]; ++vi[1])
            for (vi[2]=L[4]; vi[2]<=H[4]; ++vi[2]) {
                V = velocity(vi);
                d = Distribution(i,j,vi[0],vi[1],vi[2]);
                Jx(i,j) += V[0]*d;
                Jy(i,j) += V[1]*d;
                Jz(i,j) += V[2]*d;
                
                Vxx(i,j) += V[0]*V[0]*d;
                Vxy(i,j) += V[0]*V[1]*d;
                Vxz(i,j) += V[0]*V[2]*d;
                Vyy(i,j) += V[1]*V[1]*d;
                Vyz(i,j) += V[1]*V[2]*d;
                Vzz(i,j) += V[2]*V[2]*d;
                
            }
      }
    }
    
    boundary->ScalarFieldReduce(Jx);
    boundary->ScalarFieldReduce(Jy);
    boundary->ScalarFieldReduce(Jz);

    boundary->ScalarFieldReduce(Vxx);
    boundary->ScalarFieldReduce(Vxy);
    boundary->ScalarFieldReduce(Vxz);
    boundary->ScalarFieldReduce(Vyy);
    boundary->ScalarFieldReduce(Vyz);
    boundary->ScalarFieldReduce(Vzz);
}

template<class ForceField>
ScalarField &VlasovSpecies<ForceField>::KineticEnergy () {
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

template<class ForceField>
double VlasovSpecies<ForceField>::TotalEnergy () {
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

template<class ForceField>
VelocityD VlasovSpecies<ForceField>::getJ(int i, int j) {
    return VelocityD(Jx(i,j), Jy(i,j), Jz(i,j));
}

template<class ForceField>
FixedArray<double,6> VlasovSpecies<ForceField>::getVVTens(int i, int j) {    
    FixedArray<double,6> Result;
    
    Result[0] = Vxx(i,j);
    Result[1] = Vxy(i,j);
    Result[2] = Vxz(i,j);
    Result[3] = Vyy(i,j);
    Result[4] = Vyz(i,j);
    Result[5] = Vzz(i,j);
    
    return Result;
    
}





template<class ForceField>
void VlasovSpecies<ForceField>::advance(double timestep) {
  switch (RKState) {
      case -2: RKState = -1;
               break;
      case -1: RKState = 0;
              advanceStepA(timestep); 
              break;
      case 0: advanceStepFull(timestep); 
              break;
  }
//    const int *UBound = Distribution.getHigh();
//    const int *LBound = Distribution.getLow();
// 
//    switch (RKState) {
//        case -2: RKState = -1;
//                 break;
//        case -1: RKState = 0;
//                for (int i=LBound[0]; i<=UBound[0]; ++i)
//                  for (int j=LBound[1]; j<=UBound[1]; ++j)
//                    for (int k=LBound[2]; k<=UBound[2]; ++k)
//                      for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                        for (int m=LBound[4]; m<=UBound[4]; ++m)
//                          TempDist(i,j,k,l,m) = Distribution(i,j,k,l,m);
//                advanceStepA(timestep/3.); 
//                break;
//        case 0: advanceStepB(timestep/3.);         // Distribution = c0 = d1
// 
//                for (int i=LBound[0]; i<=UBound[0]; ++i)
//                  for (int j=LBound[1]; j<=UBound[1]; ++j)
//                    for (int k=LBound[2]; k<=UBound[2]; ++k)
//                      for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                        for (int m=LBound[4]; m<=UBound[4]; ++m)
//                          T1(i,j,k,l,m) = Distribution(i,j,k,l,m);
// 
// 
//                advanceStepA(2*timestep/3.); 
//                RKState = 1;
//                break;
//        case 1: advanceStepB(2*timestep/3.);        // Distribution = c1
// 
//                for (int i=LBound[0]; i<=UBound[0]; ++i)
//                  for (int j=LBound[1]; j<=UBound[1]; ++j)
//                    for (int k=LBound[2]; k<=UBound[2]; ++k)
//                      for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                        for (int m=LBound[4]; m<=UBound[4]; ++m) {
//                          Distribution(i,j,k,l,m) 
//                            = Distribution(i,j,k,l,m) + TempDist(i,j,k,l,m) - T1(i,j,k,l,m);
//                          T2(i,j,k,l,m) = Distribution(i,j,k,l,m); // = d2
//                        }
//                advanceStepA(3*timestep/4.);
// 
//                RKState = 2;
//                break;
//        case 2: advanceStepB(3*timestep/4.);      // Distribution = c2
//                for (int i=LBound[0]; i<=UBound[0]; ++i)
//                  for (int j=LBound[1]; j<=UBound[1]; ++j)
//                    for (int k=LBound[2]; k<=UBound[2]; ++k)
//                      for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                        for (int m=LBound[4]; m<=UBound[4]; ++m) {
//                          Distribution(i,j,k,l,m) 
//                              =  (1/4.)*TempDist(i,j,k,l,m) 
//                               + (3/4.)*T1(i,j,k,l,m) 
//                               +  Distribution(i,j,k,l,m)
//                               - T2(i,j,k,l,m); 
//                          TempDist(i,j,k,l,m) = Distribution(i,j,k,l,m);
//                        } 
//                advanceStepA(timestep/3.); 
//                RKState = 0;
//                break;
//    }
}




template<class ForceField>
void VlasovSpecies<ForceField>::InterpolationInitStep(const VlasovDist &Dist) {
    f_infty = 1;
}

template<class ForceField>
double VlasovSpecies<ForceField>::densityError() {
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

template<class ForceField>
void VlasovSpecies<ForceField>::correctDensityError(double err) {
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    for (int i=LBound[0]; i<=UBound[0]; ++i)
      for (int j=LBound[1]; j<=UBound[1]; ++j)
        for (int k=LBound[2]; k<=UBound[2]; ++k) 
          for (int l=LBound[3]; l<=UBound[3]; ++l) 
            for (int m=LBound[4]; m<=UBound[4]; ++m) 
              Distribution(i,j,k,l,m) /= err;
}



template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepA(double timestep) {
        
    InterpolationInitStep(Distribution);
    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepB(double timestep) {

    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_z(timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);

    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
    
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepFull(double timestep) {

    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_z(timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);

    advanceSpace_x(timestep);
    advanceSpace_y(timestep);
    
}

/** @brief Calls the advance method of the Advancer base class.
 * then calls Task::Execute to execute sub-Tasks
 */
template<class ForceField>
void VlasovSpecies<ForceField>::Execute () {
    if ( (t%20) == 0) {
        double err = densityError();
        correctDensityError(err);
    }
    t++;

//    cerr << "Potential("<<t<<")\n";
//    pPot->Execute(dt);
    
//    cerr << "Advancing("<<t<<")\n";
    advance(dt);
//    cerr << "Done Advancing("<<t<<")\n";

//    write_Distribution(Distribution,"Advanced.out");
}

