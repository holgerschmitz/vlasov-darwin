#include "wrapvlasov.h"
#include <fstream>
#include "darwin.h"

// -------------------------------------------------------------------
// EFieldForce
// -------------------------------------------------------------------

void EFieldForce::Init(double dttx_, double vnorm, Potential *pPot_) {
    pE = pPot_;
    dttx = dttx_;
}

VelocityD EFieldForce::Force(const PositionI &Pos, const VelocityD &Vel) {
    double Fx = pE->GetEx()(Pos[0],Pos[1]);
    double Fy = pE->GetEy()(Pos[0],Pos[1]);
    return VelocityD(dttx*Fx, dttx*Fy, 0.0);
        
}


// -------------------------------------------------------------------
// EBFieldForce
// -------------------------------------------------------------------

void EBFieldForce::Init(double dttx_, double vnorm_, Potential *pPot_) {
    pPot = pPot_;
    dttx = dttx_;
    vnorm = vnorm_;
}

void EBFieldForce::setBField(VelocityD B_) { B = B_; }

ScalarField &EBFieldForce::GetEx() { return pPot->GetEx(); }
ScalarField &EBFieldForce::GetEy() { return pPot->GetEy(); }


VelocityD EBFieldForce::Force(const PositionI &Pos, const VelocityD &Vel) {
    double Ex = pPot->GetEx()(Pos[0],Pos[1]);
    double Ey = pPot->GetEy()(Pos[0],Pos[1]);
    const double Ez = 0;
    
    double BMag2 = B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
    double BMag = sqrt(BMag2);
    double ebx = B[0]/BMag;
    double eby = B[1]/BMag;
    double ebz = B[2]/BMag;
    
    double Epar = ebx*Ex+eby*Ey+ebz*Ez;
    double Eparx = Epar*ebx;
    double Epary = Epar*eby;
    double Eparz = Epar*ebz;

    // normalizing velocity
    double vx = Vel[0]/vnorm;
    double vy = Vel[1]/vnorm;
    double vz = Vel[2]/vnorm;

    double vpar = ebx*vx+eby*vy+ebz*vz;
    double vparx = vpar*ebx;
    double vpary = vpar*eby;
    double vparz = vpar*ebz;
    
    double vdrx = Ey*B[2]-Ez*B[1];
    double vdry = Ez*B[0]-Ex*B[2];
    double vdrz = Ex*B[1]-Ey*B[0];
    
    double vperpx = vx - vparx - vdrx;
    double vperpy = vy - vpary - vdry;
    double vperpz = vz - vparz - vdrz;
    
    double vperp = sqrt(vperpx*vperpx + vperpy*vperpy + vperpz*vperpz);
    
    double eperpx = vperpx/vperp;
    double eperpy = vperpy/vperp;
    double eperpz = vperpz/vperp;
    
    double enormx = eby*eperpz - ebz*eperpy;
    double enormy = ebz*eperpx - ebx*eperpz;
    double enormz = ebx*eperpy - eby*eperpx;
    
    double Fx = Eparx - vperp*( 0.5*eperpx*BMag2*dttx + enormx*BMag );
    double Fy = Epary - vperp*( 0.5*eperpy*BMag2*dttx + enormy*BMag );
    double Fz = Eparz - vperp*( 0.5*eperpz*BMag2*dttx + enormz*BMag );
    
    return VelocityD(dttx*Fx, dttx*Fy, dttx*Fz);
        
}


// -------------------------------------------------------------------
// EMDarwinForce
// -------------------------------------------------------------------

void EMDarwinForce::Init(double dttx_, double vnorm_, Darwin* pFields_) {
    pFields = pFields_;
    dttx = dttx_;
    vnorm = vnorm_;
}

ScalarField &EMDarwinForce::GetEx() { return pFields->GetEx(); }
ScalarField &EMDarwinForce::GetEy() { return pFields->GetEy(); }
ScalarField &EMDarwinForce::GetEz() { return pFields->GetEz(); }
ScalarField &EMDarwinForce::GetBx() { return pFields->GetBx(); }
ScalarField &EMDarwinForce::GetBy() { return pFields->GetBy(); }
ScalarField &EMDarwinForce::GetBz() { return pFields->GetBz(); }


VelocityD EMDarwinForce::Force(const PositionI &Pos, const VelocityD &Vel) {
    double Ex = pFields->GetEx()(Pos[0],Pos[1]);
    double Ey = pFields->GetEy()(Pos[0],Pos[1]);
    double Ez = pFields->GetEz()(Pos[0],Pos[1]);

    VelocityD B;
    B[0] = pFields->GetBx()(Pos[0],Pos[1]);
    B[1] = pFields->GetBy()(Pos[0],Pos[1]);
    B[2] = pFields->GetBz()(Pos[0],Pos[1]);
    
    double BMag2 = B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
    double BMag = sqrt(BMag2);
    double ebx = B[0]/BMag;
    double eby = B[1]/BMag;
    double ebz = B[2]/BMag;
    
    double Epar = ebx*Ex+eby*Ey+ebz*Ez;
    double Eparx = Epar*ebx;
    double Epary = Epar*eby;
    double Eparz = Epar*ebz;

    // normalizing velocity
    double vx = Vel[0]/vnorm;
    double vy = Vel[1]/vnorm;
    double vz = Vel[2]/vnorm;

    double vpar = ebx*vx+eby*vy+ebz*vz;
    double vparx = vpar*ebx;
    double vpary = vpar*eby;
    double vparz = vpar*ebz;
    
    double vdrx = Ey*B[2]-Ez*B[1];
    double vdry = Ez*B[0]-Ex*B[2];
    double vdrz = Ex*B[1]-Ey*B[0];
    
    double vperpx = vx - vparx - vdrx;
    double vperpy = vy - vpary - vdry;
    double vperpz = vz - vparz - vdrz;
    
    double vperp = sqrt(vperpx*vperpx + vperpy*vperpy + vperpz*vperpz);
    
    double eperpx = vperpx/vperp;
    double eperpy = vperpy/vperp;
    double eperpz = vperpz/vperp;
    
    double enormx = eby*eperpz - ebz*eperpy;
    double enormy = ebz*eperpx - ebx*eperpz;
    double enormz = ebx*eperpy - eby*eperpx;
    
    double Fx = Eparx - vperp*( 0.5*eperpx*BMag2*dttx + enormx*BMag );
    double Fy = Epary - vperp*( 0.5*eperpy*BMag2*dttx + enormy*BMag );
    double Fz = Eparz - vperp*( 0.5*eperpz*BMag2*dttx + enormz*BMag );
    
    return VelocityD(dttx*Fx, dttx*Fy, dttx*Fz);
        
}


void write_Distribution(VlasovDist &Dist, const char *fname) {
    ofstream O(fname);
    for (int i=0; i<=Dist.getHigh(1); ++i) 
      for (int j=0; j<=Dist.getHigh(0); ++j) 
        for (int k=0; k<=Dist.getHigh(2); ++k) 
          for (int l=0; l<=Dist.getHigh(3); ++l) 
            for (int m=0; m<=Dist.getHigh(4); ++m) 
                O   << i << " "  << j << " "  << k << " "  
                    << l << " "  << m << " "
                    << Dist(j,i,k,l,m) << endl;
}



