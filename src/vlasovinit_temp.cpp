// -*- C++ -*-
// $Id$

//=========================================================================
//===============   VlasovMaxwellInit  ====================================
//=========================================================================

template<class ForceField>
VlasovMaxwellInit<ForceField>::VlasovMaxwellInit(VlasovMaxwellInitData<ForceField> &data) {
  pVlasov = data.vlasov;
  
  u_stream[0] = data.InitStream_vx;
  u_stream[1] = data.InitStream_vy;
  u_stream[2] = data.InitStream_vz;
  
  v_th[0] = data.InitTherm_vx;
  v_th[1] = data.InitTherm_vy;
  v_th[2] = data.InitTherm_vz;
  
  N = data.Init_N0;
  n_perturb = data.Init_N1;
  
  k_perturb[0] = data.Init_kx;
  k_perturb[1] = data.Init_ky;
  
  cout << "VlasovMaxwellInit::Init (" << u_stream << ") (" <<  v_th << ") ";
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
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  double Nx = GlDistHigh[0]-GlDistLow[0]-3;
  double Ny = GlDistHigh[1]-GlDistLow[1]-3;
  
  PositionI Xi;
  VelocityI Vi;
  VelocityD UStream(u_stream[0],u_stream[1],u_stream[2]);
  VelocityD VTh(v_th[0],v_th[1],v_th[2]);
  
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
      
      double NPert = (
		      1+n_perturb*cos(
				      2*PI*k_perturb[0]*Xi[0]/Nx
				      +2*PI*k_perturb[1]*Xi[1]/Ny
				      )
		      );
      
      //        double NPert =  1 + n_perturb*frand();
      
      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
	for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
	  for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
	    
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
	      if (Vi[j]==L[j+2]) {
		F1[j] = 0.5*(erf(vd1p[j]) + 1);
	      } else if (Vi[j]==H[j+2]) {
		F1[j] = 0.5*(1 - erf(vd1m[j]));
	      } else {
		F1[j] = 0.5*(erf(vd1p[j]) - erf(vd1m[j]));

	      }
	    }
                
                
	    double F = N*NPert*F1.product();
	    dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) = F;
	  }
        
      // Set all the boundaries to zero
        
      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
	for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) {
	  dist(Xi[0],Xi[1],Vi[0],Vi[1],L[4])=0; 
	  dist(Xi[0],Xi[1],Vi[0],Vi[1],H[4])=0; 
        }
           
      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
	for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
	  dist(Xi[0],Xi[1],Vi[0],L[3],Vi[2])=0; 
	  dist(Xi[0],Xi[1],Vi[0],H[3],Vi[2])=0; 
        }
      
      for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
	for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
	  dist(Xi[0],Xi[1],L[2],Vi[1],Vi[2])=0; 
	  dist(Xi[0],Xi[1],H[2],Vi[1],Vi[2])=0; 
        }
    }
  //    std::cout << "Initialized " << Xi[0] << std::endl;
}


//=========================================================================
//===============   VlasovTwoMaxwellInit  =================================
//=========================================================================

template<class ForceField>
VlasovTwoMaxwellInit<ForceField>::VlasovTwoMaxwellInit(VlasovTwoMaxwellInitData<ForceField> &data) {
  pVlasov = data.vlasov;
    
  u_stream1[0] = data.InitStream_vx;
  u_stream1[1] = data.InitStream_vy;
  u_stream1[2] = data.InitStream_vz;
    
  u_stream2[0] = data.InitStream2_vx;
  u_stream2[1] = data.InitStream2_vy;
  u_stream2[2] = data.InitStream2_vz;
    
  v_th1[0] = data.InitTherm_vx;
  v_th1[1] = data.InitTherm_vy;
  v_th1[2] = data.InitTherm_vz;
    
  v_th2[0] = data.InitTherm2_vx;
  v_th2[1] = data.InitTherm2_vy;
  v_th2[2] = data.InitTherm2_vz;
    
  N1 = data.Init_N0;
  N2 = data.Init2_N0;
  n_perturb1 = data.Init_N1;
  n_perturb2 = data.Init2_N1;
    
  k_perturb1[0] = data.Init_kx;
  k_perturb1[1] = data.Init_ky;

  k_perturb2[0] = data.Init2_kx;
  k_perturb2[1] = data.Init2_ky;

  cout << "VlasovTwoMaxwellInit::Init (" <<
    u_stream1 << ") (" <<
    v_th1 << ") "<<
    u_stream2 << ") (" <<
    v_th2 << ")\n";
}

template<class ForceField>
VlasovTwoMaxwellInit<ForceField>::~VlasovTwoMaxwellInit() {}


/** @brief Do the initialisation
 *  
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
template<class ForceField>
void VlasovTwoMaxwellInit<ForceField>::initialise(VlasovDist &dist, const VelocityD &VelRange) {
  const int *L = dist.getLow();
  const int *H = dist.getHigh();

  double Nx = GlDistHigh[0]-GlDistLow[0]-3;
  double Ny = GlDistHigh[1]-GlDistLow[1]-3;

  PositionI Xi;
  VelocityI Vi;
  VelocityD UStream1(u_stream1[0],u_stream1[1],u_stream1[2]);
  VelocityD VTh1(v_th1[0],v_th1[1],v_th1[2]);
  VelocityD UStream2(u_stream2[0],u_stream2[1],u_stream2[2]);
  VelocityD VTh2(v_th2[0],v_th2[1],v_th2[2]);
        
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0]) 
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
      cerr << "X pos " << Xi << endl;
      double NPert1 = (
		       1+n_perturb1*cos(
					2*PI*k_perturb1[0]*Xi[0]/Nx
					+2*PI*k_perturb1[1]*Xi[1]/Ny
					)
		       );
        
        
      double NPert2 = (
		       1+n_perturb2*cos(
					2*PI*k_perturb2[0]*Xi[0]/Nx
					+2*PI*k_perturb2[1]*Xi[1]/Ny
					)
		       );
        

      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) {
	for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) { 
	  for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
	    //                cerr << "V pos " << Vi << endl;

	    VelocityD V( pVlasov->velocity(Vi) );
	    VelocityD Vm( 
			 (pVlasov->velocity(Vi-VelocityI(1,1,1)) + V)*0.5 
			 );
	    VelocityD Vp( 
			 (pVlasov->velocity(Vi+VelocityI(1,1,1)) + V)*0.5
			 );
                
                
	    VelocityD vd1m((Vm - UStream1)/VTh1);
	    VelocityD vd1p((Vp - UStream1)/VTh1);

	    VelocityD vd2m((Vm - UStream2)/VTh2);
	    VelocityD vd2p((Vp - UStream2)/VTh2);

	    VelocityD F1, F2;
                
	    for (int j=0; j<3; ++j) {
	      if (Vi[j]==L[j+2]) {
		F1[j] = 0.5*(erf(vd1p[j]) + 1);
		F2[j] = 0.5*(erf(vd2p[j]) + 1);
	      } else if (Vi[j]==H[j+2]) {
		F1[j] = 0.5*(1 - erf(vd1m[j]));
		F2[j] = 0.5*(1 - erf(vd2m[j]));
	      } else {
		F1[j] = 0.5*(erf(vd1p[j]) - erf(vd1m[j]));
		F2[j] = 0.5*(erf(vd2p[j]) - erf(vd2m[j]));
	      }
	    }
                
                
	    double F = N1*NPert1*F1.product() + N2*NPert2*F2.product();
	    dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) = F;
	  }
	}
      }
        
      // Set all the boundaries to zero
        
      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) {
	for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) {
	  dist(Xi[0],Xi[1],Vi[0],Vi[1],L[4])=0; 
	  dist(Xi[0],Xi[1],Vi[0],Vi[1],H[4])=0; 
	}
      }
           
      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
	for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
	  dist(Xi[0],Xi[1],Vi[0],L[3],Vi[2])=0; 
	  dist(Xi[0],Xi[1],Vi[0],H[3],Vi[2])=0; 
	}
      
      for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
	for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
	  dist(Xi[0],Xi[1],L[2],Vi[1],Vi[2])=0; 
	  dist(Xi[0],Xi[1],H[2],Vi[1],Vi[2])=0; 
	}
    }
  //    std::cout << "Initialized " << Xi[0] << std::endl;
}


//=========================================================================
//===============   VlasovWaveGenInit   ===================================
//=========================================================================

template<class ForceField>
VlasovWaveGenInit<ForceField>::VlasovWaveGenInit(VlasovMaxwellInitData<ForceField> &data) {
  pVlasov = data.vlasov;
  
  u_stream[0] = data.InitStream_vx;
  u_stream[1] = data.InitStream_vy;
  u_stream[2] = data.InitStream_vz;
  
  v_th[0] = data.InitTherm_vx;
  v_th[1] = data.InitTherm_vy;
  v_th[2] = data.InitTherm_vz;
  
  N = data.Init_N0;
  n_perturb = data.Init_N1;
    
  cout << "VlasovWaveGenInit::Init (" << v_th << ") ";

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
  const int *L = dist.getLow();
  const int *H = dist.getHigh();

  int Nx = H[0]-L[0]-3;
  int Ny = H[1]-L[1]-3;

  PositionI Xi;
  VelocityI Vi;
  VelocityD VTh(v_th[0],v_th[1],v_th[2]);
  VelocityD UStream;
    
  double NPert;
  int state=0;
    
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0]) 
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
        
//      double PP = cos( 2*PI*(Xi[0]/double(GlGridX)));
//      double PZ = sin( 2*PI*(Xi[0]/double(GlGridX)));
      double PP=0;
      double PZ=0;
      for (int kx=1; kx<(GlGridX/2); ++kx)  {
//      for (int kx=1; kx<2; ++kx)  {
     PP += cos( 2*PI*kx*(Xi[0]/double(GlGridX)-kx/15.0));
     PZ += sin( 2*PI*kx*(Xi[0]/double(GlGridX)-kx/15.0));
     //PP += cos( 2*PI*kx*(Xi[0]/double(GlGridX)));
     //PZ += sin( 2*PI*kx*(Xi[0]/double(GlGridX)));
        
      }
        
      //        for (int ky=1; ky<Ny/2; ++kx) 
      //            PP += cos( 2*PI*ky*Xi[0]/double(Ny)) / double(ky);

      PP *= 0.01;
      PZ *= 0.01;
        
      double NPert = 1.0;
        
      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
     for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
       for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {

         VelocityD V( pVlasov->velocity(Vi) );
         VelocityD Vm( 
              (pVlasov->velocity(Vi-VelocityI(1,1,1)) + V)*0.5 
              );
         VelocityD Vp( 
              (pVlasov->velocity(Vi+VelocityI(1,1,1)) + V)*0.5
              );


         UStream[0] = u_stream[0]*PP;
         UStream[1] = u_stream[1]*PP;
         UStream[2] = u_stream[2]*PZ;

         VelocityD vd1m((Vm-UStream)/VTh);

         VelocityD vd1p((Vp-UStream)/VTh);

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
      state = (state==0)?1:0;

    }
  //    std::cout << "Initialized " << Xi[0] << std::endl;

}

// //=========================================================================
// //===============   VlasovCurrentSheetInit ================================
// //=========================================================================
// 
// template<class ForceField>
// VlasovCurrentSheetInit<ForceField>::VlasovCurrentSheetInit(VlasovSpecies<ForceField>* pVlasov_) {
//   pVlasov = pVlasov_;
//   
//   u_stream[0] = InitStream_vx;
//   u_stream[1] = InitStream_vy;
//   u_stream[2] = InitStream_vz;
//   
//   v_th[0] = InitTherm_vx;
//   v_th[1] = InitTherm_vy;
//   v_th[2] = InitTherm_vz;
//   
//   N = Init_N0;
//     
//   cout << "VlasovMaxwellInit::Init (" << u_stream << ") (" <<  v_th << ") ";
// }
// 
// template<class ForceField>
// VlasovCurrentSheetInit<ForceField>::~VlasovCurrentSheetInit() {}
// 
// 
// /** @brief Do the initialisation
//  *   
//  *  Iterate through the whole distribution function and assign the appropriate
//  *  phase space density to every point in phase space. The Phase space density 
//  *  is calculated as a Maxwellian distribution.
//  */
// template<class ForceField>
// void VlasovCurrentSheetInit<ForceField>::initialise(VlasovDist &dist, const VelocityD &VelRange) {
//   const int *L = dist.getLow();
//   const int *H = dist.getHigh();
//   
//   double Nx = GlDistHigh[0]-GlDistLow[0]-3;
//   double Ny = GlDistHigh[1]-GlDistLow[1]-3;
//   
//   int currup = (GlDistHigh[0]-GlDistLow[0]-3)/4;
//   int curdown = (3*(GlDistHigh[0]-GlDistLow[0]-3))/4;
//   
//   PositionI Xi;
//   VelocityI Vi;
//   VelocityD UStream(u_stream[0],u_stream[1],u_stream[2]);
//   VelocityD VTh(v_th[0],v_th[1],v_th[2]);
//   
//   for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
//     for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
//       
//       
//       for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
// 	    for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
// 	      for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
// 	    
// 	        VelocityD V( pVlasov->velocity(Vi) );
// 	        VelocityD Vm( 
// 			     (pVlasov->velocity(Vi-VelocityI(1,1,1)) + V)*0.5 
// 			     );
// 	        VelocityD Vp( 
// 			     (pVlasov->velocity(Vi+VelocityI(1,1,1)) + V)*0.5
// 			     );
// 
// 
// 	        VelocityD vd1m,vd1p;
//             if ( (Xi[0] == currup) ) {
//               vd1m = (Vm - UStream)/VTh;
//               vd1p = (Vp - UStream)/VTh;
//             } else if ( (Xi[0] == curdown) ) {
//               vd1m = (Vm + UStream)/VTh;
//               vd1p = (Vp + UStream)/VTh;
//             } else {
//               vd1m = Vm/VTh;
//               vd1p = Vp/VTh;
//             }
// 
// 	        VelocityD F1;
// 
// 	        for (int j=0; j<3; ++j) {
// 	          if (Vi[j]==L[j+2]) {
// 		    F1[j] = 0.5*(erf(vd1p[j]) + 1);
// 	          } else if (Vi[j]==H[j+2]) {
// 		    F1[j] = 0.5*(1 - erf(vd1m[j]));
// 	          } else {
// 		    F1[j] = 0.5*(erf(vd1p[j]) - erf(vd1m[j]));
// 
// 	          }
// 	        }
// 
// 
// 	        double F = N*F1.product();
// 	        dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) = F;
// 	  }
//         
//       // Set all the boundaries to zero
//         
// /*      for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
// 	    for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) {
// 	      dist(Xi[0],Xi[1],Vi[0],Vi[1],L[4])=0; 
// 	      dist(Xi[0],Xi[1],Vi[0],Vi[1],H[4])=0; 
//         }
//            
//       for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
// 	    for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
// 	      dist(Xi[0],Xi[1],Vi[0],L[3],Vi[2])=0; 
// 	      dist(Xi[0],Xi[1],Vi[0],H[3],Vi[2])=0; 
//         }
//       
//       for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
// 	    for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
// 	      dist(Xi[0],Xi[1],L[2],Vi[1],Vi[2])=0; 
// 	      dist(Xi[0],Xi[1],H[2],Vi[1],Vi[2])=0; 
//         }*/
//     } 
//   //    std::cout << "Initialized " << Xi[0] << std::endl;
// }
// 
