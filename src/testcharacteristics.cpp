
#include "vlasov.h"

double dt = 0.01;
double Charge = -1;
VelocityD EField(0,0,0);
VelocityD BField(0,0,-1);
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

int main()
{
  VelocityD Pos = XStart;
  VelocityD Vel = VStart;
  VelocityD dV;
  for (int i=0; i<629; ++i)
  {
    dV = Force(Vel,dt);
    Vel = Vel + dV;
    Pos = Pos + Vel*dt;
    cout << i << " " << Vel << " " << Pos << "\n";
  }
}

