#include "bernstein.h"
#include <cmath>

#include <fstream>

//=========================================================================
//===============   VlasovMaxwellInit  ====================================
//=========================================================================

VlasovBernsteinInit::VlasovBernsteinInit() {}

PARAMETERMAP* VlasovBernsteinInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N"] = WParameter(new ParameterValue<double>(&N, 1));
  (*pm)["amplitude"] = WParameter(new ParameterValue<double>(&amp, 0));
  (*pm)["kx"] = WParameter(new ParameterValue<int>(&wavevec[0], 0));
  (*pm)["ky"] = WParameter(new ParameterValue<int>(&wavevec[1], 0));
  (*pm)["distk"] = WParameter(new ParameterValue<double>(&distk, 1));
  (*pm)["omega"] = WParameter(new ParameterValue<double>(&omega, 1));
  (*pm)["accuracy"] = WParameter(new ParameterValue<int>(&accuracy, 10));
  return pm;
}

VlasovBernsteinInit::~VlasovBernsteinInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void VlasovBernsteinInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  
  PositionI Xi;
  VelocityI Vi;
  VelocityD VTh(v_th[0],v_th[1],v_th[2]);

  for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
	  for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1])
    {
      Vi[2] = 0;
      VelocityD V( pVlasov->velocity(Vi) );
	    VelocityD Vm( 
			  (pVlasov->velocity(Vi-VelocityI(1,1,1)) + V)*0.5 
			);
	    VelocityD Vp( 
			  (pVlasov->velocity(Vi+VelocityI(1,1,1)) + V)*0.5
			);

	    VelocityD vd1m(Vm/VTh);

	    VelocityD vd1p(Vp/VTh);

	    Complex dist1 = Distribution1Integral(Vm[0], Vm[1], Vp[0], Vp[1], accuracy);
      
      VelocityD F1;
                
      for (int j=0; j<2; ++j) {
        if (Vi[j]==L[j+2]) {
	        F1[j] = 0.5*(erf(vd1p[j]) + 1);
        } else if (Vi[j]==H[j+2]) {
          F1[j] = 0.5*(1 - erf(vd1m[j]));
        } else {
          F1[j] = 0.5*(erf(vd1p[j]) - erf(vd1m[j]));
        }
      }
      
      for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
        for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) 
        {
          double phase 
            = 2*PIl*wavevec[0]*Xi[0]/Nx
				     +2*PIl*wavevec[1]*Xi[1]/Ny;
             
          Complex perturb = exp(Complex(0,phase));
      
     
	        for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) 
          {
	    
	          VelocityD V( pVlasov->velocity(Vi) );
	          VelocityD Vm( 
			        (pVlasov->velocity(Vi-VelocityI(1,1,1)) + V)*0.5 
			      );
	          VelocityD Vp( 
			        (pVlasov->velocity(Vi+VelocityI(1,1,1)) + V)*0.5
			      );
	    
                
	          vd1m[2] = Vm[2]/VTh[2];

	          vd1p[2] = Vp[2]/VTh[2];
               
	          if (Vi[2]==L[4]) {
		          F1[2] = 0.5*(erf(vd1p[2]) + 1);
            } else if (Vi[2]==H[4]) {
      		    F1[2] = 0.5*(1 - erf(vd1m[2]));
      	    } else {
		          F1[2] = 0.5*(erf(vd1p[2]) - erf(vd1m[2]));
      	    }
                
	          double F = N*F1.product() + real(perturb*dist1)*F1[2];
	          dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) = F;
	        }
        }
    }
      // Set all the boundaries to zero
    for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
        for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1])
        {
          for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
	          for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
            {
	            dist(Xi[0],Xi[1],Vi[0],Vi[1],L[4])=0; 
	            dist(Xi[0],Xi[1],Vi[0],Vi[1],H[4])=0; 
            }
           
          for (Vi[0] = L[2]; Vi[0] <= H[2]; ++Vi[0]) 
	          for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) 
            {
	            dist(Xi[0],Xi[1],Vi[0],L[3],Vi[2])=0; 
	            dist(Xi[0],Xi[1],Vi[0],H[3],Vi[2])=0; 
            }

          for (Vi[1] = L[3]; Vi[1] <= H[3]; ++Vi[1]) 
            for (Vi[2] = L[4]; Vi[2] <= H[4]; ++Vi[2]) 
            {
	            dist(Xi[0],Xi[1],L[2],Vi[1],Vi[2])=0; 
	            dist(Xi[0],Xi[1],H[2],Vi[1],Vi[2])=0; 
            }
        }
  double err = pVlasov->densityError(dist)/N;
  pVlasov->correctDensityError(err, dist);
  
  //    std::cout << "Initialized " << Xi[0] << std::endl;
}

VlasovBernsteinInit::Complex 
  VlasovBernsteinInit::ComplexCore(double vperp, double phi, int n)
{
  double sp = sin(phi);
  double cp = cos(phi);
  Complex argu(0,-n*phi + distk*vperp*sp);
  Complex mult(sp, -(n-omega)*cp);
  return mult*exp(argu);
}


double VlasovBernsteinInit::DistFactor(double vperp, int n)
{
  double nw = n-omega;
  return vperp*jn(n,distk*vperp)*exp(-vperp*vperp)/(nw*nw-1);
}

VlasovBernsteinInit::Complex 
  VlasovBernsteinInit::Distribution1(double vperp, double phi)
{
  std::list<Complex> summands;
    
  double max;
  double absterm;
  int n=1;

  Complex term = DistFactor(vperp,0)*ComplexCore(vperp,phi,0);
  summands.push_back(term);
  absterm = abs(term);
  max=absterm;
  
  while (absterm/max > 1e-14)
  {
    term = DistFactor(vperp,n)*ComplexCore(vperp,phi,n) + DistFactor(vperp,-n)*ComplexCore(vperp,phi,-n);
    summands.push_back(term);
    ++n;
    absterm = abs(term);
    if (absterm>max) max=absterm;
  }
  
  Complex sum=0;

  for 
  (
    std::list<Complex>::iterator it = summands.begin(); 
    it != summands.end();
    ++it
  )
  {
    sum += *it;
  }
    
  return sum;
}

VlasovBernsteinInit::Complex 
  VlasovBernsteinInit::Distribution1Integral
  (
    double vxmin, double vymin, 
    double vxmax, double vymax,
    int steps
  )
{
  double dvx = (vxmax-vxmin)/double(steps);
  double dvy = (vymax-vymin)/double(steps);
  
  Complex sum = 0;
  
  double x = vxmin - 0.5*dvx;
  for (int i=0; i<steps; ++i)
  {
    x +=dvx;
    double y = vymin - 0.5*dvy;
    for (int j=0; j<steps; ++j)
    {
      y +=dvy;
      
      double v = sqrt(x*x+y*y);
      double phi = atan2(y,x);
      
      sum += Distribution1(v,phi);
    }
  }
  return amp*sum*dvx*dvy;
}
