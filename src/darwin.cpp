// -*- C++ -*-
// $Id$

#include "darwin.h"
#include "globals.h"
#include "process.h"

// ----------------------------------------------------------------------
// Darwin

void Darwin::AddSpecies(EMDarwinForce* pS) { 	
    species.push_back(pEMDarwinForce(pS));
}

void Darwin::Init () {
  tstep= 0;
  n0 = -1;
  double VRatio = Parameters::instance().velocityRatio();
  csc = VRatio*VRatio;
  
  //  Fields range from 0 to N+1, where the inner points are
  //  1 to N
  
  LBound = Process::instance().getBoundary().scalarLow();
  HBound = Process::instance().getBoundary().scalarHigh();
  
  dx[0] = Parameters::instance().gridSpace_x();
  dx[1] = Parameters::instance().gridSpace_y();
  
  dV = Parameters::instance().volumeQuant();
	
  std::cout << "Grid Spacing is " << dx << endl;
  std::cout << "Grid Size is " << LBound << " to " << HBound << endl;
    
  // resize grid
  den.resize(LBound.Data(),HBound.Data());
  om2.resize(LBound.Data(),HBound.Data());
    
  jx.resize(LBound.Data(),HBound.Data());
  jy.resize(LBound.Data(),HBound.Data());
  jz.resize(LBound.Data(),HBound.Data());

  sx.resize(LBound.Data(),HBound.Data());
  sy.resize(LBound.Data(),HBound.Data());
  sz.resize(LBound.Data(),HBound.Data());

  vxx.resize(LBound.Data(),HBound.Data());
  vxy.resize(LBound.Data(),HBound.Data());
  vxz.resize(LBound.Data(),HBound.Data());
  vyy.resize(LBound.Data(),HBound.Data());
  vyz.resize(LBound.Data(),HBound.Data());
  vzz.resize(LBound.Data(),HBound.Data());

  Pot.resize(LBound.Data(),HBound.Data());
  Az.resize(LBound.Data(),HBound.Data());

  Ex.resize(LBound.Data(),HBound.Data());
  Ey.resize(LBound.Data(),HBound.Data());
  Ez.resize(LBound.Data(),HBound.Data());
    
  Etx.resize(LBound.Data(),HBound.Data());
  Ety.resize(LBound.Data(),HBound.Data());
  Theta.resize(LBound.Data(),HBound.Data());
  DivEt.resize(LBound.Data(),HBound.Data());

  Bx.resize(LBound.Data(),HBound.Data());
  By.resize(LBound.Data(),HBound.Data());
  Bz.resize(LBound.Data(),HBound.Data());
    
  den.clear();
  om2.clear();
    
  jx.clear();
  jy.clear();
  jz.clear();

  sx.clear();
  sy.clear();
  sz.clear();

  Pot.clear();
  Az.clear();

  Ex.clear();
  Ey.clear();
  Ez.clear();

  Etx.clear();
  Ety.clear();
  Theta.clear();
  DivEt.clear();

  Bx.clear();
  By.clear();
  Bz.clear();

  int gssx = HBound[0]-LBound[0], gssy = HBound[1]-LBound[1];
  pois = new Poisson;
  pois->resize(
    PositionD(0.0,0.0), 
    PositionD(
      gssx*Parameters::instance().gridSpace_x(), 
      gssy*Parameters::instance().gridSpace_y()
    ), 
    PositionI(gssx,gssy),
    Poisson::periodic, 
    Poisson::periodic
  );

  helmh = new Helmholtz();

  In.resize(LBound.Data(),HBound.Data());
  Lambda.resize(LBound.Data(),HBound.Data());
  Out.resize(LBound.Data(),HBound.Data());

  std::cout << "Done Darwin: Size=( " << LBound << "),(" << HBound<< ")" << endl;
}

bool Darwin::Execute () {
    
  int i;
  double dF, dF2;
  
  int lx0 = LBound[0], lx1 = LBound[0]+1;
  int ly0 = LBound[1], ly1 = LBound[1]+1;
  int mx0 = HBound[0], mx1 = HBound[0]-1;
  int my0 = HBound[1], my1 = HBound[1]-1;
  
  ScalarField tmp;
  
  tstep++;

  tmp.resize(LBound.Data(),HBound.Data());
  
  //    std::cout << "Initialising Density\n";
  // initialise
    
  /* *************************************
   *  Clearing densities and current densities
   */
  den.clear();
  om2.clear();
    
  jx.clear();
  jy.clear();
  jz.clear();

  sx.clear();
  sy.clear();
  sz.clear();

  vxx.clear();
  vxy.clear();
  vxz.clear();
  vyy.clear();
  vyz.clear();
  vzz.clear();
	
  VelocityD jt;
  FixedArray<double,6> vvt;
    
  /* *************************************
   *  First ask the species to create the density and the current density
   *  and add them all up.
   */
    
  // iterate through the species
  for (int s = species.size() - 1; s >= 0; s--) {
    EMDarwinForce* pS = species[s];

    dF = pS->getCharge()/dV;
    dF2 = dF * pS->getCharge();
        
    //        cerr << "Creating Density\n";
    DistMomentRho *distRho = pS->getDerivedRho();
    DistMomentVelocities *distVel = pS->getDerivedVelocities();
    
    // ... and get it
    ScalarField &rho = distRho->getRho();
    for (int j=ly0; j<=my0; ++j) 
      for (int i=lx0; i<=mx0; ++i) {
	    den(i,j) += dF*rho(i,j);
	    om2(i,j) += dF2*rho(i,j);

	    jt = distVel->getJ(i,j);
	    jx(i,j) += dF*jt[0];
	    jy(i,j) += dF*jt[1];
	    jz(i,j) += dF*jt[2];

	    sx(i,j) += dF2*jt[0];
	    sy(i,j) += dF2*jt[1];
	    sz(i,j) += dF2*jt[2];

	    vvt =  distVel->getVVTens(i,j);
	    vxx(i,j) +=  dF*vvt[0];
	    vxy(i,j) +=  dF*vvt[1];
	    vxz(i,j) +=  dF*vvt[2];
	    vyy(i,j) +=  dF*vvt[3];
	    vyz(i,j) +=  dF*vvt[4];
	    vzz(i,j) +=  dF*vvt[5];
      }
  }

  /// @todo Diagnostic should be moved to separate classes
//   if(mainproc /*&& ((tstep%10)==0)*/) {
//     //if (false) {
//     write_Scalar(den, "den.out");
//     write_Scalar(om2, "om2.out");
//     write_Scalar(jx, "jx.out");
//     write_Scalar(jy, "jy.out");
//     write_Scalar(jz, "jz.out");
//     write_Scalar(sx, "sx.out");
//     write_Scalar(sy, "sy.out");
//     write_Scalar(sz, "sz.out");
//     write_Scalar(vxx, "vxx.out");
//     write_Scalar(vxy, "vxy.out");
//     write_Scalar(vxz, "vxz.out");
//     write_Scalar(vyy, "vyy.out");
//     write_Scalar(vyz, "vyz.out");
//     write_Scalar(vzz, "vzz.out");
//   }

  /* *************************************
   *  With the charge density we can first calculate the 
   *  scalar potential Pot. 
   *  By differentiating we get the longitudinal electric field.
   *  This is stored provisionally in the Ex and Ey fields.
   *  Remember: Ez = 0 in two dimensions
   */

  for (int j=ly0; j<=my0; ++j) 
    for (int i=lx0; i<=mx0; ++i) { 
      In(i,j) = -(den(i,j)+n0);
      //            std::cerr << "den " << i << " " << j << " " << In(i,j) << std::endl;
    }
    
  pois->solve(In,Pot);
    
  /// @todo Diagnostic should be moved to separate classes
//  if(mainproc && ((tstep%100)==0)) write_Scalar(Pot, "Pot.out");
  //if (false) write_Scalar(Pot, "Pot.out");

  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      Ex(i,j) = (Pot(i-1,j) - Pot(i+1,j)) / (2*dx[0]);

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      Ey(i,j) = (Pot(i,j-1) - Pot(i,j+1)) / (2*dx[1]);

  // Perform wrapping
  for (int j=ly0; j<=my0; ++j) {
    Ex(lx0,j) = Ex(mx1,j);
    Ex(mx0,j) = Ex(lx1,j);
  }

  // Perform wrapping
  for (int i=lx0; i<=mx0; ++i) {
    Ey(i,ly0) = Ey(i,my1);
    Ey(i,my0) = Ey(i,ly1);
  }


  /* *************************************
   *  Now we want to calculate the magnetic field B.
   *  In two dimensions we calculate the z-component of the
   *  vector potential A_z which leads to the magnetic field
   *  components B_x and B_y.
   *  The component B_z can be calculated directly from rot j_z.
   */

  /* *************************************
   * The z--component of the vector potential.
   */
    
  pois->solve(jz,Az);

    
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      Az(i,j) *= -csc;
            
  /* *************************************
   * ... resulting in Bx and By
   */

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      Bx(i,j) = (Az(i,j+1) - Az(i,j-1)) / (2*dx[1]);

  for (int i=lx0; i<=mx0; ++i) {
    Bx(i,ly0) = Bx(i,my1);
    Bx(i,my0) = Bx(i,ly1);
  }
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      By(i,j) = (Az(i-1,j) - Az(i+1,j)) / (2*dx[0]);

  for (int j=ly0; j<=my0; ++j) {
    By(lx0,j) = By(mx1,j);
    By(mx0,j) = By(lx1,j);
  }

  /* *************************************
   * The z--component of  -rot(j).
   */

  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      In(i,j) = csc*(- (jy(i+1,j) - jy(i-1,j)) / (2*dx[0])
		     + (jx(i,j+1) - jx(i,j-1)) / (2*dx[1]));

  for (int i=lx0; i<=mx0; ++i) {
    In(i,ly0) = In(i,my1);
    In(i,my0) = In(i,lx1);
  }

  for (int j=ly0; j<=my0; ++j) {
    In(lx0,j) = In(mx1,j);
    In(mx0,j) = In(lx1,j);
  }

  double rotjzsum = 0;
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) {
      rotjzsum += In(i,j);
  }
  rotjzsum /= double( (mx0-lx1)*(my0-ly1) );
  
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      In(i,j) -= rotjzsum;
  }

  pois->solve(In,Bz);
    
  double Bzsum = 0;
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) {
      Bzsum += Bz(i,j);
  }
  Bzsum /= double( (mx0-lx1)*(my0-ly1) );
  
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      Bz(i,j) -= Bzsum;
      //==================
      // ACHTUNG TEST
      //==================
      //Bz(i,j) = 0;
  }
  
  /* *************************************
   *  Finally we add the global magnetic fields
   *  to the solution for B
   */
 
  const VelocityD &GlB = Parameters::instance().bField();
  for (int i=lx0; i<=mx0; ++i) {
    for (int j=ly0; j<=my0; ++j) {
      Bx(i,j) += GlB[0];
      By(i,j) += GlB[1];
      Bz(i,j) += GlB[2];
    }
  }
     
  /* *************************************
   *  The transverse electric field has to be calculated individually
   *  for each component. 
   */
    
  // This is the same for every component
  for (int i=lx0; i<=mx0; ++i) {
    for (int j=ly0; j<=my0; ++j) {
      Lambda(i,j) = csc*om2(i,j);
    }
  }
           
  // x-component of transverse electric field
  for (int i=lx1; i<=mx1; ++i) {
    for (int j=ly1; j<=my1; ++j) {
      In(i,j) = csc*(
		    //  -( vxx(i+1,j) - vxx(i-1,j)
			//  +vxy(i+1,j) - vxy(i-1,j)
			//  +vxz(i+1,j) - vxz(i-1,j)
			//  ) / (2*dx[0])                       /// -grad (rho <vv>)
		     -(vxx(i+1,j) - vxx(i-1,j)) / (2*dx[0])    
			 -(vxy(i,j+1) - vxy(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
		     +om2(i,j)*Ex(i,j)                   /// om2*E
		     +sy(i,j)*Bz(i,j)-sz(i,j)*By(i,j)    /// q rho <v> x B
		     );
    }
  }
    
  for (int i=lx0; i<=mx0; ++i) {
    In(i,ly0) = In(i,my1);
    In(i,my0) = In(i,lx1);
  }

  for (int j=ly0; j<=my0; ++j) {
    In(lx0,j) = In(mx1,j);
    In(mx0,j) = In(lx1,j);
  }

  helmh->solve(Out,In,Lambda);
    
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      Etx(i,j) = -Out(i,j);
    
  // y-component of transverse electric field
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) {
      In(i,j) = csc*(
		     //-( vxy(i,j+1) - vxy(i,j-1)
			 // +vyy(i,j+1) - vyy(i,j-1)
			 // +vyz(i,j+1) - vyz(i,j-1)
			 // ) / (2*dx[1])                       /// -grad (rho <vv>)
		     -(vxy(i+1,j) - vxy(i-1,j)) / (2*dx[0])    
			 -(vyy(i,j+1) - vyy(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
		     +om2(i,j)*Ey(i,j)                   /// om2*E
		     +sz(i,j)*Bx(i,j)-sx(i,j)*Bz(i,j)    /// q rho <v> x B
		     );
    }
    
  for (int i=lx0; i<=mx0; ++i) {
    In(i,ly0) = In(i,my1);
    In(i,my0) = In(i,lx1);
  }

  for (int j=ly0; j<=my0; ++j) {
    In(lx0,j) = In(mx1,j);
    In(mx0,j) = In(lx1,j);
  }

  helmh->solve(Out,In,Lambda);
    
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      Ety(i,j) = -Out(i,j);
    }
    
  // z-component of transverse electric field
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) {
      In(i,j) = csc*(
		     -(vxz(i+1,j) - vxz(i-1,j)) / (2*dx[0])    
			 -(vyz(i,j+1) - vyz(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
             +sx(i,j)*By(i,j)-sy(i,j)*Bx(i,j));/// q rho <v> x B
    }
    
  for (int i=lx0; i<=mx0; ++i) {
    In(i,ly0) = In(i,my1);
    In(i,my0) = In(i,lx1);
  }

  for (int j=ly0; j<=my0; ++j) {
    In(lx0,j) = In(mx1,j);
    In(mx0,j) = In(lx1,j);
  }

  helmh->solve(Out,In,Lambda);
    
  //==================
  // ACHTUNG TEST
  //==================
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      Ez(i,j) = -Out(i,j);
      //Ez(i,j) = 0;
    }
           
  clearDiv(Etx, Ety);

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      Ex(i,j) += Etx(i,j);
      Ey(i,j) += Ety(i,j);
    }
        
  /* *************************************
   *  This should be it!
   *  Now we can write out the diagnostics
   */
    
  /// @todo Diagnostic should be moved to separate classes
//   if(mainproc /*&& ((tstep%100)==0)*/)   {
//     //if (false) {
//     write_Scalar(Ex, "Ex.out");
//     write_Scalar(Ey, "Ey.out");
//     write_Scalar(Ez, "Ez.out");
//     write_Scalar(Bx, "Bx.out");
//     write_Scalar(By, "By.out");
//     write_Scalar(Bz, "Bz.out");
//     write_Scalar(Az, "Az.out");
//     write_Scalar(Pot,"Phi.out");
//   } 

  //    Task::Execute();
  return false;
    
}

void Darwin::clearDiv(ScalarField &Fx, ScalarField &Fy) {
  int lx0 = LBound[0], lx1 = LBound[0]+1;
  int ly0 = LBound[1], ly1 = LBound[1]+1;
  int mx0 = HBound[0], mx1 = HBound[0]-1;
  int my0 = HBound[1], my1 = HBound[1]-1;
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j)  
      DivEt(i,j) = (Fx(i+1,j) - Fx(i-1,j) ) /(2*dx[0])
	+(Fy(i,j+1) - Fy(i,j-1) ) /(2*dx[1]);
                        
  for (int i=lx0; i<=mx0; ++i) {
    DivEt(i,ly0) = DivEt(i,my1);
    DivEt(i,my0) = DivEt(i,lx1);
  }

  for (int j=ly0; j<=my0; ++j) {
    DivEt(lx0,j) = DivEt(mx1,j);
    DivEt(mx0,j) = DivEt(lx1,j);
  }

  pois->solve(DivEt,Theta);
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) { 
      Fx(i,j) -= (Theta(i+1,j) - Theta(i-1,j) ) /(2*dx[0]);
      Fy(i,j) -= (Theta(i,j+1) - Theta(i,j-1) ) /(2*dx[1]);
    }
    
  for (int i=lx0; i<=mx0; ++i) {
    Fx(i,ly0) = Fx(i,my1);
    Fx(i,my0) = Fx(i,lx1);
    Fy(i,ly0) = Fy(i,my1);
    Fy(i,my0) = Fy(i,lx1);
  }

  for (int j=ly0; j<=my0; ++j) {
    Fx(lx0,j) = Fx(mx1,j);
    Fx(mx0,j) = Fx(lx1,j);
    Fy(lx0,j) = Fy(mx1,j);
    Fy(mx0,j) = Fy(lx1,j);
  }


}
