#include "derivedfields.h"
#include "vlasovbase.h"

DistributionDerivedField::DistributionDerivedField(Boundary *boundary_) 
    : boundary(boundary_) {}
    
void DerivedFieldsContainer::add(pDistributionDerivedField field) {
  std::string key((*field).name());
  derivedFields[key] = field;
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
  cerr << "Field Class " << name << " not found\n";
  return pDistributionDerivedField(NULL);
}

DistMomentRho::DistMomentRho(Boundary *boundary_) 
    : DistributionDerivedField(boundary_) 
{
  PositionI Lowx = boundary->scalarLow();
  PositionI Highx = boundary->scalarHigh();
  
  Rho.resize(Lowx.Data(),Highx.Data());  
}
    

void DistMomentRho::calc(ForceFieldBase &vlasov) {

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
    
    boundary->ScalarFieldReduce(Rho); 
}

ScalarField& DistMomentRho::getField(std::string name) {
  return getRho();
}

DistMomentVelocities::DistMomentVelocities(Boundary *boundary_) 
    : DistributionDerivedField(boundary_) 
{
  PositionI Lowx = boundary->scalarLow();
  PositionI Highx = boundary->scalarHigh();
  
  Jx.resize(Lowx.Data(),Highx.Data());
  Jy.resize(Lowx.Data(),Highx.Data());
  Jz.resize(Lowx.Data(),Highx.Data());

  Vxx.resize(Lowx.Data(),Highx.Data());
  Vxy.resize(Lowx.Data(),Highx.Data());
  Vxz.resize(Lowx.Data(),Highx.Data());
  Vyy.resize(Lowx.Data(),Highx.Data());
  Vyz.resize(Lowx.Data(),Highx.Data());
  Vzz.resize(Lowx.Data(),Highx.Data());
}
    

void DistMomentVelocities::calc(ForceFieldBase &vlasov) {
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

VelocityD DistMomentVelocities::getJ(int i, int j) {
    return VelocityD(Jx(i,j), Jy(i,j), Jz(i,j));
}



FixedArray<double,6> DistMomentVelocities::getVVTens(int i, int j) {    
    FixedArray<double,6> Result;
    
    Result[0] = Vxx(i,j);
    Result[1] = Vxy(i,j);
    Result[2] = Vxz(i,j);
    Result[3] = Vyy(i,j);
    Result[4] = Vyz(i,j);
    Result[5] = Vzz(i,j);
    
    return Result;
    
}

ScalarField& DistMomentVelocities::getField(std::string name) {
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
