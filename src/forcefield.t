// -------------------------------------------------------------------
// GenericEMForce
// -------------------------------------------------------------------

template<class FieldType>
void GenericEMForce<FieldType>::Init(double dttx_) {
    dttx = dttx_;
    Rho = new DistMomentRho(boundary);
    Veloc = new DistMomentVelocities(boundary);
    derivedFields.add(Rho);
    derivedFields.add(Veloc);
}

template<class FieldType>
ScalarField &GenericEMForce<FieldType>::GetEx() { return pPot->GetEx(); }

template<class FieldType>
ScalarField &GenericEMForce<FieldType>::GetEy() { return pPot->GetEy(); }

template<class FieldType>
ScalarField &GenericEMForce<FieldType>::GetEz() { return pPot->GetEz(); }

template<class FieldType>
ScalarField &GenericEMForce<FieldType>::GetBx() { return pPot->GetBx(); }

template<class FieldType>
ScalarField &GenericEMForce<FieldType>::GetBy() { return pPot->GetBy(); }

template<class FieldType>
ScalarField &GenericEMForce<FieldType>::GetBz() { return pPot->GetBz(); }


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
    double Ex = Charge*GetEx()(Pos[0],Pos[1])/Mass;
    double Ey = Charge*GetEy()(Pos[0],Pos[1])/Mass;
    double Ez = Charge*GetEz()(Pos[0],Pos[1])/Mass;

    double Bx = Charge*GetBx()(Pos[0],Pos[1])/Mass;
    double By = Charge*GetBy()(Pos[0],Pos[1])/Mass;
    double Bz = Charge*GetBz()(Pos[0],Pos[1])/Mass;

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

template<class FieldType>
ScalarField &GenericEMForce<FieldType>::FieldEnergy() {
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

template<class FieldType>
DistMomentRho *GenericEMForce<FieldType>::getDerivedRho() {
  return Rho;
}

template<class FieldType>
DistMomentVelocities *GenericEMForce<FieldType>::getDerivedVelocities() {
  return Veloc;
}


