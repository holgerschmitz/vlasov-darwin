#include "potential.h"
#include "darwin.h"

inline double GenericEMForceBase_Electrostatic::GetEx(const PositionI &Pos) 
{ 
  return pPot->GetEx()(Pos[0],Pos[1]); 
}

inline double GenericEMForceBase_Electrostatic::GetEy(const PositionI &Pos)
{ 
  return pPot->GetEy()(Pos[0],Pos[1]); 
}


// -------------------------------------------------------------------
// GenericEMForce
// -------------------------------------------------------------------

template<class FieldType>
void GenericEMForceBase_FullEM<FieldType>::Init(double dttx_) {
    dttx = dttx_;
    Rho = new DistMomentRho(boundary);
    Veloc = new DistMomentVelocities(boundary);
    derivedFields.add(Rho);
    derivedFields.add(Veloc);
}

template<class FieldType>
inline double GenericEMForceBase_FullEM<FieldType>::GetEx(const PositionI &Pos) 
{ return pPot->GetEx()(Pos[0],Pos[1]); }

template<class FieldType>
inline double GenericEMForceBase_FullEM<FieldType>::GetEy(const PositionI &Pos) 
{ return pPot->GetEy()(Pos[0],Pos[1]); }

template<class FieldType>
inline double GenericEMForceBase_FullEM<FieldType>::GetEz(const PositionI &Pos) 
{ return pPot->GetEz()(Pos[0],Pos[1]); }

template<class FieldType>
inline double GenericEMForceBase_FullEM<FieldType>::GetBx(const PositionI &Pos) 
{ return pPot->GetBx()(Pos[0],Pos[1]); }

template<class FieldType>
inline double GenericEMForceBase_FullEM<FieldType>::GetBy(const PositionI &Pos) 
{ return pPot->GetBy()(Pos[0],Pos[1]); }

template<class FieldType>
inline double GenericEMForceBase_FullEM<FieldType>::GetBz(const PositionI &Pos) 
{ return pPot->GetBz()(Pos[0],Pos[1]); }


template<class FieldType>
DistMomentRho *GenericEMForceBase_FullEM<FieldType>::getDerivedRho() {
  return Rho;
}

template<class FieldType>
DistMomentVelocities *GenericEMForceBase_FullEM<FieldType>::getDerivedVelocities() {
  return Veloc;
}

template<class FieldType>
VelocityD GenericEMForce<FieldType>
                       ::Force(const PositionI &Pos, 
                               const VelocityD &Vel,
                               double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;

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

template<class FieldType>
double GenericEMForce<FieldType>
                       ::ForceX(const PositionI &Pos, 
                                const VelocityD &Vel,
                                double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;
    
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
    double sy = sfact*ty;
    double sz = sfact*tz;
    
    // b) now v-prime
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    double vprz = Vmz + Vmx*ty-Vmy*tx;
    
    // c) and finally V-plus
    double Vpx = Vmx + vpry*sz-vprz*sy;
    
    // Calculate new velocity minus old velocity
    double Vdiffx = Vpx + 0.5*Ex*dt - vx;
    
   
    return Vdiffx;
        
}

template<class FieldType>
double GenericEMForce<FieldType>
                       ::ForceY(const PositionI &Pos, 
                                const VelocityD &Vel,
                                double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;

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
    double sz = sfact*tz;
    
    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vprz = Vmz + Vmx*ty-Vmy*tx;
    
    // c) and finally V-plus
    double Vpy = Vmy + vprz*sx-vprx*sz;
    
    // Calculate new velocity minus old velocity
    double Vdiffy = Vpy + 0.5*Ey*dt - vy;
    
   
    return Vdiffy;
        
}

template<class FieldType>
double GenericEMForce<FieldType>
                       ::ForceZ(const PositionI &Pos, 
                                const VelocityD &Vel,
                                double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;

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
    
    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    
    // c) and finally V-plus
    double Vpz = Vmz + vprx*sy-vpry*sx;
    
    // Calculate new velocity minus old velocity
    double Vdiffz = Vpz + 0.5*Ez*dt - vz;
    
   
    return Vdiffz;
        
}

//template<class FieldType>
//ScalarField &GenericEMForce<FieldType>::FieldEnergy() {
//    ScalarField &Ex = GetEx();
//    ScalarField &Ey = GetEy();
//    ScalarField &Ez = GetEz();
//    ScalarField &Bx = GetBx();
//    ScalarField &By = GetBy();
//    ScalarField &Bz = GetBz();
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
//            FEngy(i,j) += sqr(Ez(i,j));
//            FEngy(i,j) += sqr(Bx(i,j));
//            FEngy(i,j) += sqr(By(i,j));
//            FEngy(i,j) += sqr(Bz(i,j));
//    }
//    
//    return FEngy;
//}

//********************************************************************
//********** GenericEMForceBoris2  **************************
//********************************************************************

template<class FieldType>
VelocityD GenericEMForceBoris2<FieldType>
                       ::Force(const PositionI &Pos, 
                               const VelocityD &Vel,
                               double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;

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

template<class FieldType>
double GenericEMForceBoris2<FieldType>
                       ::ForceX(const PositionI &Pos, 
                                const VelocityD &Vel,
                                double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;
    
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
    double sy = sfact*ty;
    double sz = sfact*tz;
    
    // b) now v-prime
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    double vprz = Vmz + Vmx*ty-Vmy*tx;
    
    // c) and finally V-plus
    double Vpx = Vmx + vpry*sz-vprz*sy;
    
    // Calculate new velocity minus old velocity
    double Vdiffx = Vpx + 0.5*Ex*dt - vx;
    
   
    return Vdiffx;
        
}

template<class FieldType>
double GenericEMForceBoris2<FieldType>
                       ::ForceY(const PositionI &Pos, 
                                const VelocityD &Vel,
                                double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;

    // a) Calculate t and s
    double tx = 0.5*Bx*dt;
    double ty = 0.5*By*dt;
    double tz = 0.5*Bz*dt;
    
    double sfact = 2.0/(1 + tx*tx + ty*ty + tz*tz);
    double sx = sfact*tx;
    double sy = sfact*ty;
    double sz = sfact*tz;

    // Here we assume that vx is actually given at time t+1
    // so we first have to calculate Vmx from it

    // Calculate V-minus
    double Vpx = vx-0.5*Ex*dt;
    double Vmy = vy+0.5*Ey*dt;
    double Vmz = vz+0.5*Ez*dt;
    
    double Ax = 1-sy*ty-sz*tz;
    double Vmx = (Vpx - Vmy*(sz+sy*tx) + Vmz*(sy-sz*tx)) / Ax;
    
    // Rotate
    
    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vprz = Vmz + Vmx*ty-Vmy*tx;
    
    // c) and finally V-plus
    double Vpy = Vmy + vprz*sx-vprx*sz;
    
    // Calculate new velocity minus old velocity
    double Vdiffy = Vpy + 0.5*Ey*dt - vy;
    
   
    return Vdiffy;
        
}

template<class FieldType>
double GenericEMForceBoris2<FieldType>
                       ::ForceZ(const PositionI &Pos, 
                                const VelocityD &Vel,
                                double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*GetEx(Pos)/Mass;
    double Ey = Charge*GetEy(Pos)/Mass;
    double Ez = Charge*GetEz(Pos)/Mass;

    double Bx = Charge*GetBx(Pos)/Mass;
    double By = Charge*GetBy(Pos)/Mass;
    double Bz = Charge*GetBz(Pos)/Mass;

    // Rotate
    // a) Calculate t and s
    double tx = 0.5*Bx*dt;
    double ty = 0.5*By*dt;
    double tz = 0.5*Bz*dt;
    
    double sfact = 2.0/(1 + tx*tx + ty*ty + tz*tz);
    double sx = sfact*tx;
    double sy = sfact*ty;
    double sz = sfact*tz;
    
    // Here we assume that vx and vy are actually given at time t+1
    // so we first have to replace them with the old values

    // Calculate V-minus
    double Vpx = vx-0.5*Ex*dt;
    double Vpy = vy-0.5*Ey*dt;
    double Vmz = vz+0.5*Ez*dt;
    
    double ax = 1-sy*ty-sz*tz;
    double ay = 1-sx*tx-sz*tz;
    double bx = sz + sy*tx;
    double by = sz - sx*ty;
    double cx = Vpx + (sy-sz*tx)*Vmz;
    double cy = Vpy - (sx+sz*ty)*Vmz;
    double S = (ax*ay+bx*by);
    
    double Vmx = (ay*cx-bx*cy)/S;
    double Vmy = (by*cx+ax*cy)/S;    

    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    
    // c) and finally V-plus
    double Vpz = Vmz + vprx*sy-vpry*sx;
    
    // Calculate new velocity minus old velocity
    double Vdiffz = Vpz + 0.5*Ez*dt - vz;
    
   
    return Vdiffz;
        
}


//********************************************************************
//********** GenericEMForceDirect  **************************
//********************************************************************

template<class FieldType>
VelocityD GenericEMForceDirect<FieldType>
                       ::Force(const PositionI &Pos, 
                               const VelocityD &Vel,
                               double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    double af = dt*Charge/Mass;
    
    // Storing E and B field
    double Ex = GetEx(Pos);
    double Ey = GetEy(Pos);
    double Ez = GetEz(Pos);

    double Bx = GetBx(Pos);
    double By = GetBy(Pos);
    double Bz = GetBz(Pos);

    // Calculate Force
    double Fx = Ex + vy*Bz - vz*By;
    double Fy = Ey + vz*Bx - vx*Bz;
    double Fz = Ez + vx*By - vy*Bx;
       
    return VelocityD(Fx*af, Fy*af, Fz*af);
        
}
