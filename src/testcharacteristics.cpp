
#include "vlasov.h"

double dt = 0.05;
double Charge = 1;
VelocityD EField(0,0,0);
VelocityD BField(1,0,0);
VelocityD VStart(0,1,0);
VelocityD XStart(0,0,0);

VelocityD Force(const VelocityD &Vel,
                double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*EField[0];
    double Ey = Charge*EField[1];
    double Ez = Charge*EField[2];

    double Bx = Charge*BField[0];
    double By = Charge*BField[1];
    double Bz = Charge*BField[2];
    
    // Calculate V-minus
    double Vmx = vx+0.5*Ex*dt;
    double Vmy = vy+0.5*Ey*dt;
    double Vmz = vz+0.5*Ez*dt;
    
//    cout << "B1 Vm " << Vmx << " " << Vmy << " " << Vmz << "\n";
    
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

double ForceX(const VelocityD &Vel,
                 double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*EField[0];
    double Ey = Charge*EField[1];
    double Ez = Charge*EField[2];

    double Bx = Charge*BField[0];
    double By = Charge*BField[1];
    double Bz = Charge*BField[2];
    
    // Calculate V-minus
    double Vmx = vx+0.5*Ex*dt;
    double Vmy = vy+0.5*Ey*dt;
    double Vmz = vz+0.5*Ez*dt;

//    cout << "B2x Vm " << Vmx << " " << Vmy << " " << Vmz << "\n";
    
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

double ForceY(const VelocityD &Vel,
                 double dt) {
   
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*EField[0];
    double Ey = Charge*EField[1];
    double Ez = Charge*EField[2];

    double Bx = Charge*BField[0];
    double By = Charge*BField[1];
    double Bz = Charge*BField[2];

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
    double Vmx = (Vpx - Vmy*(sz+sy*tx) - Vmz*(sy-sz*tx)) / Ax;

//    cout << "B2y Vm " << Vmx << " " << Vmy << " " << Vmz << "\n";
    
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

double ForceZ(const VelocityD &Vel,
              double dt) {
    
    // normalizing velocity
    double vx = Vel[0];
    double vy = Vel[1];
    double vz = Vel[2];
    
    // Storing E and B field
    double Ex = Charge*EField[0];
    double Ey = Charge*EField[1];
    double Ez = Charge*EField[2];

    double Bx = Charge*BField[0];
    double By = Charge*BField[1];
    double Bz = Charge*BField[2];

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

//    cout << "B2z Vm " << Vmx << " " << Vmy << " " << Vmz << "\n";

    // b) now v-prime
    double vprx = Vmx + Vmy*tz-Vmz*ty;
    double vpry = Vmy + Vmz*tx-Vmx*tz;
    
    // c) and finally V-plus
    double Vpz = Vmz + vprx*sy-vpry*sx;
    
    // Calculate new velocity minus old velocity
    double Vdiffz = Vpz + 0.5*Ez*dt - vz;
    
   
    return Vdiffz;
        
}

int main()
{
  VelocityD Pos1 = XStart;
  VelocityD Vel1 = VStart;
  VelocityD Pos2 = XStart;
  VelocityD Vel2 = VStart;
  VelocityD dV;
  double d1v;



  for (int i=0; i<600; ++i)
  {
    dV = Force(Vel1,dt);
    Vel1 = Vel1 + dV;
    Pos1 = Pos1 + Vel1*dt;

    d1v = ForceX(Vel2,dt);
    Vel2[0] += d1v;
//    cout << "("  << Vel2 << ")\n";
    d1v = ForceY(Vel2,dt);
    Vel2[1] += d1v;
//    cout << "("  << Vel2 << ")\n";
    d1v = ForceZ(Vel2,dt);
    Vel2[2] += d1v;
//    cout << "("  << Vel2 << ")\n";
    Pos2 = Pos2 + Vel2*dt;

    cout << Vel2 <<  "\n";
  }
 
}

