#include "forcefield.h"
#include "potential.h"
#include "darwin.h"

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
    double Fx = GetEx()(Pos[0],Pos[1]);
    double Fy = GetEy()(Pos[0],Pos[1]);
    return VelocityD(dttx*dt*Fx, dttx*dt*Fy, 0.0);
        
}

ScalarField &EFieldForce::FieldEnergy() {
    ScalarField &Ex = GetEx();
    ScalarField &Ey = GetEy();
    
    const int *L = Ex.getLow();
    const int *H = Ex.getHigh();
    
    FEngy.resize(Ex);
    FEngy.clear();
    
    for (int i=L[0]; i<=H[0]; ++i)
        for (int j=L[1]; i<=H[1]; ++i) {
            FEngy(i,j) += sqr(Ex(i,j));
            FEngy(i,j) += sqr(Ey(i,j));
    }
    
    return FEngy;
}

DistMomentRho *EFieldForce::getDerivedRho() {
  return Rho;
}


// -------------------------------------------------------------------
// EBFieldForce
// -------------------------------------------------------------------

void EBFieldForce::Init(double dttx_) {
  dttx = dttx_;
  Rho = new DistMomentRho(boundary);
  derivedFields.add(Rho);
}

void EBFieldForce::setBField(VelocityD B_) { B = B_; }

ScalarField &EBFieldForce::GetEx() { return pPot->GetEx(); }
ScalarField &EBFieldForce::GetEy() { return pPot->GetEy(); }


VelocityD EBFieldForce::Force(const PositionI &Pos, 
                               const VelocityD &Vel,
                               double dt) {
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx()(Pos[0],Pos[1]);
    double Ey = Charge*GetEy()(Pos[0],Pos[1]);
    double Ez = 0;

    double Bx = Charge*B[0];
    double By = Charge*B[1];
    double Bz = Charge*B[2];
    
    // Calculate V-minus
    double Vmx = vx+0.5*Ex*dt;
    double Vmy = vy+0.5*Ey*dt;
    double Vmz = vz+0.5*Ez*dt;
    
    // Rotate
    // a) Calculate t and s
    double tx = 0.5*Bx*dt;
    double ty = 0.5*By*dt;
    double tz = 0.5*Bz*dt;
    
    double sfact = 2.0/(1 + tx*tx + ty*ty + tz*tz);
    double sx = sfact*tx;
    double sy = sfact*ty;
    double sz = sfact*tz;
    
    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    double vprz = Vmz + Vmx*ty-Vmy*tx;
    
    // c) and finally V-plus
    double Vpx = Vmx + vpry*sz-vprz*sy;
    double Vpy = Vmy + vprz*sx-vprx*sz;
    double Vpz = Vmz + vprx*sy-vpry*sx;
    
    // Calculate new velocity minus old velocity
    double Vdiffx = Vpx + 0.5*Ex*dt - vx;
    double Vdiffy = Vpy + 0.5*Ey*dt - vy;
    double Vdiffz = Vpz + 0.5*Ez*dt - vz;
    
   
    return VelocityD(Vdiffx, Vdiffy, Vdiffz);
    
}

DistMomentRho *EBFieldForce::getDerivedRho() {
  return Rho;
}

// -------------------------------------------------------------------
// ConstEBFieldForce
// -------------------------------------------------------------------

void ConstEBFieldForce::Init(double dttx_) {
    dttx = dttx_;
    Veloc = new DistMomentVelocities(boundary);
    derivedFields.add(Veloc);
}

void ConstEBFieldForce::setEField(VelocityD E_) { E = E_; }

void ConstEBFieldForce::setBField(VelocityD B_) { B = B_; }

VelocityD ConstEBFieldForce::Force(const PositionI &Pos, 
                                   const VelocityD &Vel,
                                   double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*E[0];
    double Ey = Charge*E[1];
    double Ez = Charge*E[2];

    double Bx = Charge*B[0];
    double By = Charge*B[1];
    double Bz = Charge*B[2];

    // Calculate V-minus
    double Vmx = vx+0.5*Ex*dt;
    double Vmy = vy+0.5*Ey*dt;
    double Vmz = vz+0.5*Ez*dt;
    
    // Rotate
    // a) Calculate t and s
    double tx = 0.5*Bx*dt;
    double ty = 0.5*By*dt;
    double tz = 0.5*Bz*dt;
    
    double sfact = 2.0/(1 + tx*tx + ty*ty + tz*tz);
    double sx = sfact*tx;
    double sy = sfact*ty;
    double sz = sfact*tz;
    
    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    double vprz = Vmz + Vmx*ty-Vmy*tx;
    
    // c) and finally V-plus
    double Vpx = Vmx + vpry*sz-vprz*sy;
    double Vpy = Vmy + vprz*sx-vprx*sz;
    double Vpz = Vmz + vprx*sy-vpry*sx;
    
    // Calculate new velocity minus old velocity
    double Vdiffx = Vpx + 0.5*Ex*dt - vx;
    double Vdiffy = Vpy + 0.5*Ey*dt - vy;
    double Vdiffz = Vpz + 0.5*Ez*dt - vz;
    
   
    return VelocityD(Vdiffx, Vdiffy, Vdiffz);
        
}

DistMomentVelocities *ConstEBFieldForce::getDerivedVelocities() {
  return Veloc;
}


// -------------------------------------------------------------------
// EMDarwinForce
// -------------------------------------------------------------------

void EMDarwinForce::Init(double dttx_) {
    dttx = dttx_;
    Rho = new DistMomentRho(boundary);
    Veloc = new DistMomentVelocities(boundary);
    derivedFields.add(Rho);
    derivedFields.add(Veloc);
}

ScalarField &EMDarwinForce::GetEx() { return pPot->GetEx(); }
ScalarField &EMDarwinForce::GetEy() { return pPot->GetEy(); }
ScalarField &EMDarwinForce::GetEz() { return pPot->GetEz(); }
ScalarField &EMDarwinForce::GetBx() { return pPot->GetBx(); }
ScalarField &EMDarwinForce::GetBy() { return pPot->GetBy(); }
ScalarField &EMDarwinForce::GetBz() { return pPot->GetBz(); }


VelocityD EMDarwinForce::Force(const PositionI &Pos, 
                               const VelocityD &Vel,
                               double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx()(Pos[0],Pos[1]);
    double Ey = Charge*GetEy()(Pos[0],Pos[1]);
    double Ez = Charge*GetEz()(Pos[0],Pos[1]);

    double Bx = Charge*GetBx()(Pos[0],Pos[1]);
    double By = Charge*GetBy()(Pos[0],Pos[1]);
    double Bz = Charge*GetBz()(Pos[0],Pos[1]);

//    Ex = 0; Ey = 0; Ez = 0;
//    Bx = 0; By = 0; Bz = 1;
    
    // Calculate V-minus
    double Vmx = vx+0.5*Ex*dt;
    double Vmy = vy+0.5*Ey*dt;
    double Vmz = vz+0.5*Ez*dt;
    
    // Rotate
    // a) Calculate t and s
    double tx = 0.5*Bx*dt;
    double ty = 0.5*By*dt;
    double tz = 0.5*Bz*dt;
    
    double sfact = 2.0/(1 + tx*tx + ty*ty + tz*tz);
    double sx = sfact*tx;
    double sy = sfact*ty;
    double sz = sfact*tz;
    
    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    double vprz = Vmz + Vmx*ty-Vmy*tx;
    
    // c) and finally V-plus
    double Vpx = Vmx + vpry*sz-vprz*sy;
    double Vpy = Vmy + vprz*sx-vprx*sz;
    double Vpz = Vmz + vprx*sy-vpry*sx;
    
    // Calculate new velocity minus old velocity
    double Vdiffx = Vpx + 0.5*Ex*dt - vx;
    double Vdiffy = Vpy + 0.5*Ey*dt - vy;
    double Vdiffz = Vpz + 0.5*Ez*dt - vz;
    
   
    return VelocityD(Vdiffx, Vdiffy, Vdiffz);
        
}

ScalarField &EMDarwinForce::FieldEnergy() {
    ScalarField &Ex = GetEx();
    ScalarField &Ey = GetEy();
    ScalarField &Ez = GetEz();
    ScalarField &Bx = GetBx();
    ScalarField &By = GetBy();
    ScalarField &Bz = GetBz();
    
    const int *L = Ex.getLow();
    const int *H = Ex.getHigh();
    
    FEngy.resize(Ex);
    FEngy.clear();
    
    for (int i=L[0]; i<=H[0]; ++i)
        for (int j=L[1]; i<=H[1]; ++i) {
            FEngy(i,j) += sqr(Ex(i,j));
            FEngy(i,j) += sqr(Ey(i,j));
            FEngy(i,j) += sqr(Ez(i,j));
            FEngy(i,j) += sqr(Bx(i,j));
            FEngy(i,j) += sqr(By(i,j));
            FEngy(i,j) += sqr(Bz(i,j));
    }
    
    return FEngy;
}

DistMomentRho *EMDarwinForce::getDerivedRho() {
  return Rho;
}

DistMomentVelocities *EMDarwinForce::getDerivedVelocities() {
  return Veloc;
}


void write_Distribution(VlasovDist &Dist, const char *fname) {
    ofstream O(fname);
    for (int i=Dist.getLow(0); i<=Dist.getHigh(0); ++i) 
      for (int j=Dist.getLow(1); j<=Dist.getHigh(1); ++j) 
        for (int k=Dist.getLow(2); k<=Dist.getHigh(2); ++k) 
          for (int l=Dist.getLow(3); l<=Dist.getHigh(3); ++l) 
            for (int m=Dist.getLow(4); m<=Dist.getHigh(4); ++m) 
                O   << i << " "  << j << " "  << k << " "  
                    << l << " "  << m << " "
                    << Dist(i,j,k,l,m) << endl;
}


