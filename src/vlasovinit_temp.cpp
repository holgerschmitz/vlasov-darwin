
//=========================================================================
//===============   VlasovMaxwellInit   ===================================
//=========================================================================

template<class ForceField>
VlasovMaxwellInit<ForceField>::VlasovMaxwellInit(VlasovSpecies<ForceField>* pVlasov_) {
    pVlasov = pVlasov_;
    
    u_stream[0] = InitStream_vx;
    u_stream[1] = InitStream_vy;
    u_stream[2] = InitStream_vz;
    
    v_th[0] = InitTherm_vx;
    v_th[1] = InitTherm_vy;
    v_th[2] = InitTherm_vz;
    
    N = Init_N0;
    n_perturb = Init_N1;
    
    k_perturb[0] = Init_kx;
    k_perturb[1] = Init_ky;

    Nx[0] = GlGridSize_x;
    Nx[1] = GlGridSize_y;


    cout << "VlasovMaxwellInit::Init (" <<
        u_stream << ") (" <<
        v_th << ") ";
}

template<class ForceField>
VlasovMaxwellInit<ForceField>::~VlasovMaxwellInit() {}


/** @brief Do the initialisation
 *  
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
template<class ForceField>
void VlasovMaxwellInit<ForceField>::initialise(VlasovDist &dist, const VelocityD &VelRange) {
    PositionI Xi;
    VelocityI Vi;
    VelocityD UStream(u_stream[0],u_stream[1],u_stream[2]);
    VelocityD VTh(v_th[0],v_th[1],v_th[2]);
        
    for (Xi[0] = 0; Xi[0] <= dist.getHigh(0); ++Xi[0]) 
      for (Xi[1] = 0; Xi[1] <= dist.getHigh(1); ++Xi[1]) {

        double NPert = (
            1+n_perturb*cos(
                2*PI*k_perturb[0]*Xi[0]/(dist.getHigh(0)-1.)
               +2*PI*k_perturb[1]*Xi[1]/(dist.getHigh(1)-1.)
            )
        );
        for (Vi[0] = 0; Vi[0] <= dist.getHigh(2); ++Vi[0]) 
          for (Vi[1] = 0; Vi[1] <= dist.getHigh(3); ++Vi[1]) 
            for (Vi[2] = 0; Vi[2] <= dist.getHigh(4); ++Vi[2]) {

                VelocityD V( pVlasov->velocity(Vi) );
                VelocityD Vm( 
                    (pVlasov->velocity(Vi-VelocityI(1,1,1)) + V)*0.5 
                );
                VelocityD Vp( 
                    (pVlasov->velocity(Vi+VelocityI(1,1,1)) + V)*0.5
                );
                
                
                VelocityD vd1m((Vm - UStream)/VTh);

                VelocityD vd1p((Vp - UStream)/VTh);

                VelocityD F1;
                
                for (int j=0; j<3; ++j) {
                    if (Vi[j]==0) {
                        F1[j] = 0.5*(erf(vd1p[j]) + 1);
                    } else if (Vi[j]==dist.getHigh(j+2)) {
                        F1[j] = 0.5*(1 - erf(vd1m[j]));
                    } else {
                        F1[j] = 0.5*(erf(vd1p[j]) - erf(vd1m[j]));

                    }
                }
                
                
                double F = N*NPert*F1.product();
                dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) = F;
          }
    }
//    std::cout << "Initialized " << Xi[0] << std::endl;

}


//=========================================================================
//===============   VlasovWaveGenInit   ===================================
//=========================================================================

template<class ForceField>
VlasovWaveGenInit<ForceField>::VlasovWaveGenInit(VlasovSpecies<ForceField>* pVlasov_) {
    pVlasov = pVlasov_;
    
    v_th[0] = InitTherm_vx;
    v_th[1] = InitTherm_vy;
    v_th[2] = InitTherm_vz;
    
    N = Init_N0;
    n_perturb = Init_N1;
    
    Nx[0] = GlGridSize_x;
    Nx[1] = GlGridSize_y;


    cout << "VlasovWaveGenInit::Init (" <<
        v_th << ") ";

}

template<class ForceField>
VlasovWaveGenInit<ForceField>::~VlasovWaveGenInit() {}


/** @brief Do the initialisation
 *  
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
template<class ForceField>
void VlasovWaveGenInit<ForceField>::initialise(VlasovDist &dist, const VelocityD &VelRange) {
    PositionI Xi;
    VelocityI Vi;
    VelocityD VTh(v_th[0],v_th[1],v_th[2]);
    
    double xmh  = dist.getHigh(0)/2.;
    double xmhh = dist.getHigh(0)/4.;
    double ymh  = dist.getHigh(1)/2.;
    double ymhh = dist.getHigh(1)/4.;
    
    double NPert;
    
    for (Xi[0] = 0; Xi[0] <= dist.getHigh(0); ++Xi[0]) 
      for (Xi[1] = 0; Xi[1] <= dist.getHigh(1); ++Xi[1]) {
        
        NPert = 0;
        if ( (Xi[0]>=8) && (Xi[0]<16) && (Xi[1]>=8) && (Xi[1]<16) )
            NPert = n_perturb;
        if ( (Xi[0]>=16) && (Xi[0]<24) && (Xi[1]>=8) && (Xi[1]<16) )
            NPert = -n_perturb;
        if ( (Xi[0]>=8) && (Xi[0]<16) && (Xi[1]>=16) && (Xi[1]<24) )
            NPert = -n_perturb;
        if ( (Xi[0]>=16) && (Xi[0]<24) && (Xi[1]>=16) && (Xi[1]<24) )
            NPert = n_perturb;

        NPert = 1+NPert;
        
        for (Vi[0] = 0; Vi[0] <= dist.getHigh(2); ++Vi[0]) 
          for (Vi[1] = 0; Vi[1] <= dist.getHigh(3); ++Vi[1]) 
            for (Vi[2] = 0; Vi[2] <= dist.getHigh(4); ++Vi[2]) {

                VelocityD V( pVlasov->velocity(Vi) );
                VelocityD Vm( 
                    (pVlasov->velocity(Vi-VelocityI(1,1,1)) + V)*0.5 
                );
                VelocityD Vp( 
                    (pVlasov->velocity(Vi+VelocityI(1,1,1)) + V)*0.5
                );
                
                
                VelocityD vd1m(Vm/VTh);

                VelocityD vd1p(Vp/VTh);

                VelocityD F1;
                
                for (int j=0; j<3; ++j) {
                    if (Vi[j]==0) {
                        F1[j] = 0.5*(erf(vd1p[j]) + 1);
                    } else if (Vi[j]==dist.getHigh(j+2)) {
                        F1[j] = 0.5*(1 - erf(vd1m[j]));
                    } else {
                        F1[j] = 0.5*(erf(vd1p[j]) - erf(vd1m[j]));

                    }
                }
                
                
                double F = N*NPert*F1.product();
                dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) = F;
          }
    }
//    std::cout << "Initialized " << Xi[0] << std::endl;

}
