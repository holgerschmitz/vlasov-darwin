#include "forcefield.h"
#include "potential.h"
#include "darwin.h"
#include "globals.h"

// -------------------------------------------------------------------
// EFieldForce
// -------------------------------------------------------------------

void EFieldForce::Init(double dttx_) {
  dttx = dttx_;
  Rho = new DistMomentRho(boundary);
  derivedFields.add(Rho);
}

ScalarField &EFieldForce::GetEx() { return pPot->GetEx(); }
ScalarField &EFieldForce::GetEy() { return pPot->GetEy(); }

VelocityD EFieldForce::Force(const PositionI &Pos, 
                             const VelocityD &Vel,
                             double dt) {
    double Fx = Charge*GetEx()(Pos[0],Pos[1])/Mass;
    double Fy = Charge*GetEy()(Pos[0],Pos[1])/Mass;
    return VelocityD(dttx*dt*Fx, dttx*dt*Fy, 0.0);
        
}

double EFieldForce::ForceX(const PositionI &Pos, 
                           const VelocityD &Vel,
                           double dt) {
    double Fx = Charge*GetEx()(Pos[0],Pos[1])/Mass;
    return dttx*dt*Fx;        
}

double EFieldForce::ForceY(const PositionI &Pos, 
                           const VelocityD &Vel,
                           double dt) {
    double Fy = Charge*GetEy()(Pos[0],Pos[1])/Mass;
    return dttx*dt*Fy;
        
}

double EFieldForce::ForceZ(const PositionI &Pos, 
                           const VelocityD &Vel,
                           double dt) {
    return 0.0;
}

//ScalarField &EFieldForce::FieldEnergy() {
//    ScalarField &Ex = GetEx();
//    ScalarField &Ey = GetEy();
//    
//    const int *L = Ex.getLow();
//    const int *H = Ex.getHigh();
//    
//    FEngy.resize(Ex);
//    FEngy.clear();
//    
//    for (int i=L[0]; i<=H[0]; ++i)
//        for (int j=L[1]; i<=H[1]; ++i) {
//            FEngy(i,j) += sqr(Ex(i,j));
//            FEngy(i,j) += sqr(Ey(i,j));
//    }
//    
//    return FEngy;
//}

DistMomentRho *EFieldForce::getDerivedRho() {
  return Rho;
}


// -------------------------------------------------------------------
// EBFieldForce
// -------------------------------------------------------------------

void GenericEMForceBase_Electrostatic::Init(double dttx_) {
  dttx = dttx_;
  Rho = new DistMomentRho(boundary);
  B = Parameters::instance().bField();
  derivedFields.add(Rho);
}


DistMomentRho *GenericEMForceBase_Electrostatic::getDerivedRho() {
  return Rho;
}

// -------------------------------------------------------------------
// ConstEBFieldForce
// -------------------------------------------------------------------

void GenericEMForceBase_ConstEB::Init(double dttx_) {
    dttx = dttx_;
    Veloc = new DistMomentVelocities(boundary);
    E = Parameters::instance().eField();
    B = Parameters::instance().bField();
    derivedFields.add(Veloc);
}



DistMomentVelocities *GenericEMForceBase_ConstEB::getDerivedVelocities() {
  return Veloc;
}



