// -*- C++ -*-
// $Id$

#include "vlasovinit.h"

#include <sstream>
#include <math.h>
#include "globals.h"
#include "hdfstream.h"

//=========================================================================
//===============   VlasovMaxwellInit  ====================================
//=========================================================================

VlasovMaxwellInit::VlasovMaxwellInit() {}

PARAMETERMAP* VlasovMaxwellInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["Stream_vx"] = WParameter(new ParameterValue<double>(&u_stream[0], 0));
  (*pm)["Stream_vy"] = WParameter(new ParameterValue<double>(&u_stream[1], 0));
  (*pm)["Stream_vz"] = WParameter(new ParameterValue<double>(&u_stream[2], 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N, 1));
  (*pm)["N1"] = WParameter(new ParameterValue<double>(&n_perturb, 0));
  (*pm)["kx"] = WParameter(new ParameterValue<int>(&k_perturb[0], 0));
  (*pm)["ky"] = WParameter(new ParameterValue<int>(&k_perturb[1], 0));
  return pm;
}

VlasovMaxwellInit::~VlasovMaxwellInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void VlasovMaxwellInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  
  PositionI Xi;
  VelocityI Vi;
  VelocityD UStream(u_stream[0],u_stream[1],u_stream[2]);
  VelocityD VTh(v_th[0],v_th[1],v_th[2]);
  
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
      
      double NPert = (
		      1+n_perturb*cos(
				      2*PIl*k_perturb[0]*Xi[0]/Nx
				      +2*PIl*k_perturb[1]*Xi[1]/Ny
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

VlasovTwoMaxwellInit::VlasovTwoMaxwellInit() {}

PARAMETERMAP* VlasovTwoMaxwellInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["Stream_vxa"] = WParameter(new ParameterValue<double>(&u_stream1[0], 0));
  (*pm)["Stream_vya"] = WParameter(new ParameterValue<double>(&u_stream1[1], 0));
  (*pm)["Stream_vza"] = WParameter(new ParameterValue<double>(&u_stream1[2], 0));
  (*pm)["Therm_vxa"] = WParameter(new ParameterValue<double>(&v_th1[0], 1));
  (*pm)["Therm_vya"] = WParameter(new ParameterValue<double>(&v_th1[1], 1));
  (*pm)["Therm_vza"] = WParameter(new ParameterValue<double>(&v_th1[2], 1));
  (*pm)["N0a"] = WParameter(new ParameterValue<double>(&N1, 1));
  (*pm)["N1a"] = WParameter(new ParameterValue<double>(&n_perturb1, 0));
  (*pm)["kxa"] = WParameter(new ParameterValue<int>(&k_perturb1[0], 0));
  (*pm)["kya"] = WParameter(new ParameterValue<int>(&k_perturb1[1], 0));
  (*pm)["Stream_vxb"] = WParameter(new ParameterValue<double>(&u_stream2[0], 0));
  (*pm)["Stream_vyb"] = WParameter(new ParameterValue<double>(&u_stream2[1], 0));
  (*pm)["Stream_vzb"] = WParameter(new ParameterValue<double>(&u_stream2[2], 0));
  (*pm)["Therm_vxb"] = WParameter(new ParameterValue<double>(&v_th2[0], 1));
  (*pm)["Therm_vyb"] = WParameter(new ParameterValue<double>(&v_th2[1], 1));
  (*pm)["Therm_vzb"] = WParameter(new ParameterValue<double>(&v_th2[2], 1));
  (*pm)["N0b"] = WParameter(new ParameterValue<double>(&N2, 1));
  (*pm)["N1b"] = WParameter(new ParameterValue<double>(&n_perturb2, 0));
  (*pm)["kxb"] = WParameter(new ParameterValue<int>(&k_perturb2[0], 0));
  (*pm)["kyb"] = WParameter(new ParameterValue<int>(&k_perturb2[1], 0));
  return pm;
}


VlasovTwoMaxwellInit::~VlasovTwoMaxwellInit() {}


/** @brief Do the initialisation
 *  
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void VlasovTwoMaxwellInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();

  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;

  PositionI Xi;
  VelocityI Vi;
  VelocityD UStream1(u_stream1[0],u_stream1[1],u_stream1[2]);
  VelocityD VTh1(v_th1[0],v_th1[1],v_th1[2]);
  VelocityD UStream2(u_stream2[0],u_stream2[1],u_stream2[2]);
  VelocityD VTh2(v_th2[0],v_th2[1],v_th2[2]);
        
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0]) 
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
//      cerr << "X pos " << Xi << endl;
      double NPert1 = (
		       1+n_perturb1*cos(
					2*PIl*k_perturb1[0]*Xi[0]/Nx
					+2*PIl*k_perturb1[1]*Xi[1]/Ny
					)
		       );
        
        
      double NPert2 = (
		      1+n_perturb2*cos(
					 2*PIl*k_perturb2[0]*Xi[0]/Nx
					 +2*PIl*k_perturb2[1]*Xi[1]/Ny
          )
		       );
        

     for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) {
	     for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) { 
	       for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) {
	       //  cerr << "V pos " << Vi << endl;

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

VlasovWaveGenInit::VlasovWaveGenInit() { }

PARAMETERMAP* VlasovWaveGenInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["Stream_vx"] = WParameter(new ParameterValue<double>(&u_stream[0], 0));
  (*pm)["Stream_vy"] = WParameter(new ParameterValue<double>(&u_stream[1], 0));
  (*pm)["Stream_vz"] = WParameter(new ParameterValue<double>(&u_stream[2], 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N, 1));
  (*pm)["N1"] = WParameter(new ParameterValue<double>(&n_perturb, 0));
  return pm;
}

VlasovWaveGenInit::~VlasovWaveGenInit() {}


/** @brief Do the initialisation
 *  
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void VlasovWaveGenInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
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
        
//      double PP = cos( 2*PIl*(Xi[0]/double(gx)));
//      double PZ = sin( 2*PIl*(Xi[0]/double(gx)));
      double PP=0;
      double PZ=0;
      double gx = Parameters::instance().gridX();
      for (int kx=1; kx<(gx/2); ++kx)  {
//      for (int kx=1; kx<2; ++kx)  {
        PP += cos( 2*PIl*kx*(Xi[0]/gx-kx/15.0));
        PZ += sin( 2*PIl*kx*(Xi[0]/gx-kx/15.0));
        //PP += cos( 2*PIl*kx*(Xi[0]/gx));
        //PZ += sin( 2*PIl*kx*(Xi[0]/gx));
        
      }
        
      //        for (int ky=1; ky<Ny/2; ++kx) 
      //            PP += cos( 2*PIl*ky*Xi[0]/double(Ny)) / double(ky);

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
//=========================================================================
//===============   VlasovGaussTempInit  ====================================
//=========================================================================

VlasovGaussTempInit::VlasovGaussTempInit() { }

PARAMETERMAP* VlasovGaussTempInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["Stream_vx"] = WParameter(new ParameterValue<double>(&u_stream[0], 0));
  (*pm)["Stream_vy"] = WParameter(new ParameterValue<double>(&u_stream[1], 0));
  (*pm)["Stream_vz"] = WParameter(new ParameterValue<double>(&u_stream[2], 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N, 1));
  (*pm)["Amplitude"] = WParameter(new ParameterValue<double>(&t_amp, 1));
  (*pm)["Width"] = WParameter(new ParameterValue<double>(&t_width, 10));
  return pm;
}

VlasovGaussTempInit::~VlasovGaussTempInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void VlasovGaussTempInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();

  int Nx = H[0]-L[0]-3;
  int Ny = H[1]-L[1]-3;
  
  double NxH = Parameters::instance().gridX() / 2;
  
  PositionI Xi;
  VelocityI Vi;
  VelocityD UStream(u_stream[0],u_stream[1],u_stream[2]);
  VelocityD VTh(v_th[0],v_th[1],v_th[2]);
  
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
      
      
      double TPert =
          1+t_amp*exp(
            -sqr( (double(Xi[0])-NxH)/t_width ) 
          );
//      std::cerr << "temperatur " << Xi[0] << " " << TPert << "\n";
      VTh[0] = sqrt(TPert)*v_th[0];
            
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


            double F = N*F1.product();
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
//========================   VlasovHDFInit ================================
//=========================================================================

VlasovHDFInit::VlasovHDFInit() { }

PARAMETERMAP* VlasovHDFInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["file"] = WParameter(new ParameterValue<std::string>(&fname, ""));
  return pm;
}

VlasovHDFInit::~VlasovHDFInit() {}

void VlasovHDFInit::initialise(ForceFieldBase *pVlasov)
{
  std::cerr << "Restart " << fname << std::endl;
  
  Parameters::instance().setRestart(true);
  
  std::string parsed=fname;

  std::ostringstream comrankstr;
  comrankstr << Process::instance().getBoundary().getUniqueId();
  std::string comrank = comrankstr.str();
 
  parsed.replace(parsed.find("#p"),2,comrank);
  
  VlasovDist &dist = pVlasov->getDistribution();  

  std::cerr << "reading " << parsed << std::endl;
  HDFistream input(parsed.c_str());
  input >> dist;
  input.close();
}


// //=========================================================================
// //===============   VlasovCurrentSheetInit ================================
// //=========================================================================
// 
// VlasovCurrentSheetInit::VlasovCurrentSheetInit(VlasovSpecies<ForceField>* pVlasov_) {
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
// void VlasovCurrentSheetInit<ForceField>::initialise(VlasovDist &dist) {
//   const int *L = dist.getLow();
//   const int *H = dist.getHigh();
//   
//  PhasePositionI GlLow  = Parameters::instance().distLow();
//  PhasePositionI GlHigh = Parameters::instance().distHigh();
//  
//  double Nx = GlHigh[0]-GlLow[0]-3;
//  double Ny = GlHigh[1]-GlLow[1]-3;
//   
//   int currup = (GlHigh[0]-GlLow[0]-3)/4;
//   int curdown = (3*(GlHigh[0]-GlLow[0]-3))/4;
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
