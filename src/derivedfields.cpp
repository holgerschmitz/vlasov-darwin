#include "derivedfields.h"
#include "vlasovbase.h"

#include <sstream>
#include "process.h"

DistributionDerivedField::DistributionDerivedField(Boundary *boundary_) 
    : boundary(boundary_) {}

DerivedFieldsContainer::DerivedFieldsContainer(Boundary *boundary_)
    : boundary(boundary_) {}

pDistributionDerivedField DerivedFieldsContainer::add(pDistributionDerivedField field) {
  std::string key((*field).name());
  cerr << "Adding field class " << key << "\n";
  derivedFields[key] = field;
  return field;
}

void DerivedFieldsContainer::update(ForceFieldBase &vlasov) {
  for (
      MapType::iterator it = derivedFields.begin();
      it != derivedFields.end();
      ++it
  )
  {
    (*(*it).second).calc(vlasov);
  }
}

pDistributionDerivedField DerivedFieldsContainer::getField(std::string name) {
  MapType::iterator pos = derivedFields.find(name);
  if (pos != derivedFields.end()) return pos->second;
  
  if ("rho"==name) {
    return add(new DistMomentRho(boundary));
  }
  if ("vels"==name) {
    return add(new DistMomentVelocities(boundary));
  }
  if ("hflux"==name) {
    return add(new DistMomentHeatFlux(boundary));
  }
  cerr << "Field Class " << name << " not found\n";
  return pDistributionDerivedField(NULL);
}

DistMomentRhoBase::DistMomentRhoBase(Boundary *boundary_) 
    : DistributionDerivedField(boundary_) 
{
  PositionI Lowx = boundary->scalarLow();
  PositionI Highx = boundary->scalarHigh();
  
  Rho.resize(Lowx.Data(),Highx.Data());
  Rho.setComponent(ScalarField::ScalarComponent);
  Rho.setParity(ScalarField::EvenParity);
}
    
ScalarField& DistMomentRhoBase::getField(std::string name) {
  return getRho();
}

void DistMomentRhoOne::calc(ForceFieldBase &vlasov) {

    VlasovDist &dist = vlasov.getDistribution();
  
    const int *L = dist.getLow();
    const int *H = dist.getHigh();
    
    Rho.clear();
    
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
                Rho(ix,iy) += dist(ix,iy,j,k,l);
//        cerr << "den " << ix << " " << iy << " " << gRho(ix,iy) << endl;
      }
    }
    
    boundary->ScalarFieldCombine(Rho); 
}

DistMomentVelocitiesBase::DistMomentVelocitiesBase(Boundary *boundary_) 
    : DistributionDerivedField(boundary_) 
{
  PositionI Lowx = boundary->scalarLow();
  PositionI Highx = boundary->scalarHigh();
  
  Jx.resize(Lowx.Data(),Highx.Data());
  Jx.setComponent(ScalarField::XComponent);
  Jx.setParity(ScalarField::OddParity);
  
  Jy.resize(Lowx.Data(),Highx.Data());
  Jy.setComponent(ScalarField::YComponent);
  Jy.setParity(ScalarField::OddParity);
  
  Jz.resize(Lowx.Data(),Highx.Data());
  Jz.setComponent(ScalarField::ZComponent);
  Jz.setParity(ScalarField::OddParity);

  Vxx.resize(Lowx.Data(),Highx.Data());
  Vxx.setComponent(ScalarField::ScalarComponent);
  Vxx.setParity(ScalarField::EvenParity);
  Vxy.resize(Lowx.Data(),Highx.Data());
  Vxy.setComponent(ScalarField::ScalarComponent);
  Vxy.setParity(ScalarField::EvenParity);
  Vxz.resize(Lowx.Data(),Highx.Data());
  Vxz.setComponent(ScalarField::ScalarComponent);
  Vxz.setParity(ScalarField::EvenParity);
  Vyy.resize(Lowx.Data(),Highx.Data());
  Vyy.setComponent(ScalarField::ScalarComponent);
  Vyy.setParity(ScalarField::EvenParity);
  Vyz.resize(Lowx.Data(),Highx.Data());
  Vyz.setComponent(ScalarField::ScalarComponent);
  Vyz.setParity(ScalarField::EvenParity);
  Vzz.resize(Lowx.Data(),Highx.Data());
  Vzz.setComponent(ScalarField::ScalarComponent);
  Vzz.setParity(ScalarField::EvenParity);
}
    


VelocityD DistMomentVelocitiesBase::getJ(int i, int j) {
    return VelocityD(Jx(i,j), Jy(i,j), Jz(i,j));
}



FixedArray<double,6> DistMomentVelocitiesBase::getVVTens(int i, int j) {    
    FixedArray<double,6> Result;
    
    Result[0] = Vxx(i,j);
    Result[1] = Vxy(i,j);
    Result[2] = Vxz(i,j);
    Result[3] = Vyy(i,j);
    Result[4] = Vyz(i,j);
    Result[5] = Vzz(i,j);
    
    return Result;
    
}

ScalarField& DistMomentVelocitiesBase::getField(std::string name) {
  if ("Jx"==name) return Jx;
  else
  if ("Jy"==name) return Jy;
  else
  if ("Jz"==name) return Jz;
  else
  if ("Vxx"==name) return Vxx;
  else
  if ("Vxy"==name) return Vxy;
  else
  if ("Vxz"==name) return Vxz;
  else
  if ("Vyy"==name) return Vyy;
  else
  if ("Vyz"==name) return Vyz;
  else
  return Vzz;
}

void DistMomentVelocitiesOne::calc(ForceFieldBase &vlasov) {
  VlasovDist &dist = vlasov.getDistribution();

  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
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
              V = vlasov.velocity(vi);
              d = dist(i,j,vi[0],vi[1],vi[2]);

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
    
  boundary->ScalarFieldCombine(Jx);
  boundary->ScalarFieldCombine(Jy);
  boundary->ScalarFieldCombine(Jz);

  boundary->ScalarFieldCombine(Vxx);
  boundary->ScalarFieldCombine(Vxy);
  boundary->ScalarFieldCombine(Vxz);
  boundary->ScalarFieldCombine(Vyy);
  boundary->ScalarFieldCombine(Vyz);
  boundary->ScalarFieldCombine(Vzz);
  
//  static int cnt=1;
//  ostringstream fname;
//  fname << "Gz"<<cnt++<<".out";
//  ofstream GzStream(fname.str().c_str());
//
//  const Boundary& boundary = Process::instance().getBoundary();
//  
//  const PositionI &LBound = boundary.scalarLow();
//  const PositionI &HBound = boundary.scalarHigh();
//
//  for (int i=LBound[0]; i<=HBound[0]; ++i) 
//    for (int j=LBound[1]; j<=HBound[1]; ++j) 
//      GzStream << i << " " << j << " " << Jz(i,j) << "\n";
//      
//  GzStream.close();
}


void DistMomentVelocitiesTwo::calc(ForceFieldBase &vlasov) {
    VlasovDist &dist = vlasov.getDistribution();

    const int *L = dist.getLow();
    const int *H = dist.getHigh();
    
    VelocityI vi;
    VelocityD V, Vm, Vp;
    double d, dpx, dmx, dpy, dmy, dpz, dmz;
    double dpxpy, dpxmy, dmxpy, dmxmy;
    double dpxpz, dpxmz, dmxpz, dmxmz;
    double dpypz, dpymz, dmypz, dmymz;
    
    Jx.clear();
    Jy.clear();
    Jz.clear();

    Vxx.clear();
    Vxy.clear();
    Vxz.clear();
    Vyy.clear();
    Vyz.clear();
    Vzz.clear();

    double dvx = vlasov.deltaVx();
    double dvy = vlasov.deltaVy();
    double dvz = vlasov.deltaVz();
    
    for (int i=L[0]+2; i<=H[0]-2; ++i) {
      for (int j=L[1]+2; j<=H[1]-2; ++j) {
        
        for (vi[0]=L[2]; vi[0]<=H[2]; ++vi[0])
          for (vi[1]=L[3]; vi[1]<=H[3]; ++vi[1])
            for (vi[2]=L[4]; vi[2]<=H[4]; ++vi[2]) {
              V = vlasov.velocity(vi);
              Vm = (vlasov.velocity(vi-VelocityI(1,1,1)) + V)*0.5;
              Vp = (vlasov.velocity(vi+VelocityI(1,1,1)) + V)*0.5;
              d = dist(i,j,vi[0],vi[1],vi[2]);
              dpx = (vi[0]<H[2])?dist(i,j,vi[0]+1,vi[1],vi[2]):0;
              dmx = (vi[0]>L[2])?dist(i,j,vi[0]-1,vi[1],vi[2]):0;
              dpy = (vi[1]<H[3])?dist(i,j,vi[0],vi[1]+1,vi[2]):0;
              dmy = (vi[1]>L[3])?dist(i,j,vi[0],vi[1]-1,vi[2]):0;
              dpz = (vi[2]<H[4])?dist(i,j,vi[0],vi[1],vi[2]+1):0;
              dmz = (vi[2]>L[4])?dist(i,j,vi[0],vi[1],vi[2]-1):0;

              dpxpy = ((vi[0]<H[2])&&(vi[1]<H[3]))?dist(i,j,vi[0]+1,vi[1]+1,vi[2]):0;
              dmxpy = ((vi[0]>L[2])&&(vi[1]<H[3]))?dist(i,j,vi[0]-1,vi[1]+1,vi[2]):0;
              dpxmy = ((vi[0]<H[2])&&(vi[1]>L[3]))?dist(i,j,vi[0]+1,vi[1]-1,vi[2]):0;
              dmxmy = ((vi[0]>L[2])&&(vi[1]>L[3]))?dist(i,j,vi[0]-1,vi[1]-1,vi[2]):0;

              dpxpz = ((vi[0]<H[2])&&(vi[2]<H[4]))?dist(i,j,vi[0]+1,vi[1],vi[2]+1):0;
              dmxpz = ((vi[0]>L[2])&&(vi[2]<H[4]))?dist(i,j,vi[0]-1,vi[1],vi[2]+1):0;
              dpxmz = ((vi[0]<H[2])&&(vi[2]>L[4]))?dist(i,j,vi[0]+1,vi[1],vi[2]-1):0;
              dmxmz = ((vi[0]>L[2])&&(vi[2]>L[4]))?dist(i,j,vi[0]-1,vi[1],vi[2]-1):0;

              dpypz = ((vi[1]<H[3])&&(vi[2]<H[4]))?dist(i,j,vi[0],vi[1]+1,vi[2]+1):0;
              dmypz = ((vi[1]>L[3])&&(vi[2]<H[4]))?dist(i,j,vi[0],vi[1]-1,vi[2]+1):0;
              dpymz = ((vi[1]<H[3])&&(vi[2]>L[4]))?dist(i,j,vi[0],vi[1]+1,vi[2]-1):0;
              dmymz = ((vi[1]>L[3])&&(vi[2]>L[4]))?dist(i,j,vi[0],vi[1]-1,vi[2]-1):0;


              Jx(i,j) += 0.25*(Vm[0]*(d+dmx) + Vp[0]*(d+dpx) );
              Jy(i,j) += 0.25*(Vm[1]*(d+dmy) + Vp[1]*(d+dpy) );
              Jz(i,j) += 0.25*(Vm[2]*(d+dmz) + Vp[2]*(d+dpz) );

              Vxx(i,j) += 0.25*(Vm[0]*Vm[0]*(d+dmx) + Vp[0]*Vp[0]*(d+dpx) );
              Vxy(i,j) += 0.0625*(    Vp[1]*Vp[0]*(d+dpx+dpy+dpxpy)
                                     +Vp[1]*Vm[0]*(d+dmx+dpy+dmxpy)
                                     +Vm[1]*Vp[0]*(d+dpx+dmy+dpxmy)
                                     +Vm[1]*Vm[0]*(d+dmx+dmy+dmxmy)
                                 );
              Vxz(i,j) += 0.0625*(    Vp[2]*Vp[0]*(d+dpx+dpz+dpxpz)
                                     +Vp[2]*Vm[0]*(d+dmx+dpz+dmxpz)
                                     +Vm[2]*Vp[0]*(d+dpx+dmz+dpxmz)
                                     +Vm[2]*Vm[0]*(d+dmx+dmz+dmxmz)
                                 );
              Vyy(i,j) += 0.25*(Vm[1]*Vm[1]*(d+dmy) + Vp[1]*Vp[1]*(d+dpy) );
              Vyz(i,j) += 0.0625*(    Vp[2]*Vp[1]*(d+dpy+dpz+dpypz)
                                     +Vp[2]*Vm[1]*(d+dmy+dpz+dmypz)
                                     +Vm[2]*Vp[1]*(d+dpy+dmz+dpymz)
                                     +Vm[2]*Vm[1]*(d+dmy+dmz+dmymz)
                                 );
              Vzz(i,j) += 0.25*(Vm[2]*Vm[2]*(d+dmz) + Vp[2]*Vp[2]*(d+dpz) );
            }
        }
    }
    
    boundary->ScalarFieldCombine(Jx);
    boundary->ScalarFieldCombine(Jy);
    boundary->ScalarFieldCombine(Jz);

    boundary->ScalarFieldCombine(Vxx);
    boundary->ScalarFieldCombine(Vxy);
    boundary->ScalarFieldCombine(Vxz);
    boundary->ScalarFieldCombine(Vyy);
    boundary->ScalarFieldCombine(Vyz);
    boundary->ScalarFieldCombine(Vzz);
}

DistMomentHeatFluxBase::DistMomentHeatFluxBase(Boundary *boundary_) 
    : DistributionDerivedField(boundary_) 
{
  PositionI Lowx = boundary->scalarLow();
  PositionI Highx = boundary->scalarHigh();
  
  HFluxX.resize(Lowx.Data(),Highx.Data());  
  HFluxX.setComponent(ScalarField::XComponent);
  HFluxX.setParity(ScalarField::OddParity);
}
    
ScalarField& DistMomentHeatFluxBase::getField(std::string name) {
  return getFluxX();
}

void DistMomentHeatFluxOne::calc(ForceFieldBase &vlasov) {
    VlasovDist &dist = vlasov.getDistribution();
  
    const int *L = dist.getLow();
    const int *H = dist.getHigh();

    VelocityI vi;
    double Vx;
    double d;
    HFluxX.clear();
        
    for (int i=L[0]+2; i<=H[0]-2; ++i)
      for (int j=L[1]+2; j<=H[1]-2; ++j) {
        
        double Jx = 0;
        for (vi[0]=L[2]; vi[0]<=H[2]; ++vi[0])
          for (vi[1]=L[3]; vi[1]<=H[3]; ++vi[1])
            for (vi[2]=L[4]; vi[2]<=H[4]; ++vi[2]) {
              Vx = vlasov.velocity(vi)[0];
              d = dist(i,j,vi[0],vi[1],vi[2]);
              Jx += Vx*d;
            }
            
        for (vi[0]=L[2]; vi[0]<=H[2]; ++vi[0])
          for (vi[1]=L[3]; vi[1]<=H[3]; ++vi[1])
            for (vi[2]=L[4]; vi[2]<=H[4]; ++vi[2]) {
              Vx = vlasov.velocity(vi)[0] - Jx;
              d = dist(i,j,vi[0],vi[1],vi[2]);
              HFluxX(i,j) += Vx*Vx*Vx*d;
            }
      }
    
    boundary->ScalarFieldCombine(HFluxX); 
}


