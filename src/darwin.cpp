// -*- C++ -*-
// $Id$

#include "darwin.h"
#include "globals.h"
#include "process.h"

#include <fstream>
// ----------------------------------------------------------------------
// Darwin

void Darwin::AddSpecies(EMDarwinForce* pS) { 	
    species.push_back(pEMDarwinForce(pS));
}

void Darwin::Init () 
{
  tstep= 0;
  n0 = Parameters::instance().bgDensity();
  double VRatio = Parameters::instance().velocityRatio();
  csc = VRatio*VRatio;
  
  //  Fields range from 0 to N+1, where the inner points are
  //  1 to N
  
  const Boundary& boundary = Process::instance().getBoundary();
  
  LBound = boundary.scalarLow();
  HBound = boundary.scalarHigh();
  
  dx[0] = Parameters::instance().gridSpace_x();
  dx[1] = Parameters::instance().gridSpace_y();
  
  dV = Parameters::instance().volumeQuant();
	
  std::cout << "Grid Spacing is " << dx << endl;
  std::cout << "Grid Size is " << LBound << " to " << HBound << endl;
    
  // resize grid
  den.resize(LBound.Data(),HBound.Data());
  den.setComponent(ScalarField::ScalarComponent);
  den.setParity(ScalarField::EvenParity);

  om2.resize(LBound.Data(),HBound.Data());
  om2.setComponent(ScalarField::ScalarComponent);
  om2.setParity(ScalarField::EvenParity);
    
  jx.resize(LBound.Data(),HBound.Data());
  jx.setComponent(ScalarField::XComponent);
  jx.setParity(ScalarField::OddParity);
  jy.resize(LBound.Data(),HBound.Data());
  jy.setComponent(ScalarField::YComponent);
  jy.setParity(ScalarField::OddParity);
  jz.resize(LBound.Data(),HBound.Data());
  jz.setComponent(ScalarField::ZComponent);
  jz.setParity(ScalarField::OddParity);

  sx.resize(LBound.Data(),HBound.Data());
  sx.setComponent(ScalarField::XComponent);
  sx.setParity(ScalarField::OddParity);
  sy.resize(LBound.Data(),HBound.Data());
  sy.setComponent(ScalarField::YComponent);
  sy.setParity(ScalarField::OddParity);
  sz.resize(LBound.Data(),HBound.Data());
  sz.setComponent(ScalarField::ZComponent);
  sz.setParity(ScalarField::OddParity);

  vxx.resize(LBound.Data(),HBound.Data());
  vxx.setComponent(ScalarField::ScalarComponent);
  vxx.setParity(ScalarField::EvenParity);

  vxy.resize(LBound.Data(),HBound.Data());
  vxy.setComponent(ScalarField::ScalarComponent);
  vxy.setParity(ScalarField::EvenParity);

  vxz.resize(LBound.Data(),HBound.Data());
  vxz.setComponent(ScalarField::ScalarComponent);
  vxz.setParity(ScalarField::EvenParity);

  vyy.resize(LBound.Data(),HBound.Data());
  vyy.setComponent(ScalarField::ScalarComponent);
  vyy.setParity(ScalarField::EvenParity);

  vyz.resize(LBound.Data(),HBound.Data());
  vyz.setComponent(ScalarField::ScalarComponent);
  vyz.setParity(ScalarField::EvenParity);

  vzz.resize(LBound.Data(),HBound.Data());
  vzz.setComponent(ScalarField::ScalarComponent);
  vzz.setParity(ScalarField::EvenParity);

  Pot.resize(LBound.Data(),HBound.Data());
  Pot.setComponent(ScalarField::ScalarComponent);
  Pot.setParity(ScalarField::EvenParity);
  
  Az.resize(LBound.Data(),HBound.Data());
  Az.setComponent(ScalarField::ZComponent);
  Az.setParity(ScalarField::OddParity);

  Ex.resize(LBound.Data(),HBound.Data());
  Ex.setComponent(ScalarField::XComponent);
  Ex.setParity(ScalarField::OddParity);
  Ey.resize(LBound.Data(),HBound.Data());
  Ey.setComponent(ScalarField::YComponent);
  Ey.setParity(ScalarField::OddParity);
  Ez.resize(LBound.Data(),HBound.Data());
  Ez.setComponent(ScalarField::ZComponent);
  Ez.setParity(ScalarField::OddParity);

    
  Etx.resize(LBound.Data(),HBound.Data());
  Etx.setComponent(ScalarField::XComponent);
  Etx.setParity(ScalarField::OddParity);
  Ety.resize(LBound.Data(),HBound.Data());
  Ety.setComponent(ScalarField::YComponent);
  Ety.setParity(ScalarField::OddParity);
  
  Theta.resize(LBound.Data(),HBound.Data());
  Theta.setComponent(ScalarField::ScalarComponent);
  Theta.setParity(ScalarField::EvenParity);
  
  DivEt.resize(LBound.Data(),HBound.Data());
  DivEt.setComponent(ScalarField::ScalarComponent);
  DivEt.setParity(ScalarField::EvenParity);

  Bx.resize(LBound.Data(),HBound.Data());
  Bx.setComponent(ScalarField::XComponent);
  Bx.setParity(ScalarField::EvenParity);
  By.resize(LBound.Data(),HBound.Data());
  By.setComponent(ScalarField::YComponent);
  By.setParity(ScalarField::EvenParity);
  Bz.resize(LBound.Data(),HBound.Data());
  Bz.setComponent(ScalarField::ZComponent);
  Bz.setParity(ScalarField::EvenParity);
    
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

  pois = new Poisson();

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
  
  const Boundary &bound = Process::instance().getBoundary();

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
//  std::cerr << "dV="<<dV<<"\n";
  for (int s = species.size() - 1; s >= 0; s--) {
    EMDarwinForce* pS = species[s];

    dF = pS->getCharge();
//    std::cerr << "dF="<<dF<<"\n";
    dF2 = dF * pS->getCharge()/pS->getMass();

    // cerr << "Creating Density\n";
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

  /* *************************************
   *  With the charge density we can first calculate the 
   *  scalar potential Pot. 
   *  By differentiating we get the longitudinal electric field.
   *  This is stored provisionally in the Ex and Ey fields.
   *  Remember: Ez = 0 in two dimensions
   */

  for (int j=ly0; j<=my0; ++j) 
    for (int i=lx0; i<=mx0; ++i) { 
      In(i,j) = (den(i,j)+n0)/csc;
      //            std::cerr << "den " << i << " " << j << " " << In(i,j) << std::endl;
    }
    
  pois->solve(Pot,In, bound.getNumBoundary(Pot));
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      Ex(i,j) = (Pot(i-1,j) - Pot(i+1,j)) / (2*dx[0]);

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      Ey(i,j) = (Pot(i,j-1) - Pot(i,j+1)) / (2*dx[1]);

  bound.ScalarFieldReduce(Ex);
  bound.ScalarFieldReduce(Ey);


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
//  std::ofstream JzStream("JZ.out");
//  for (int i=lx0; i<=mx0; ++i) 
//    for (int j=ly0; j<=my0; ++j)
//      JzStream << i << " " << j << " " << jz(i,j) << "\n";
//  JzStream.close();
  
  pois->solve(Az,jz, bound.getNumBoundary(Az));

//  std::ofstream AzStream("AZ.out");
//  for (int i=lx0; i<=mx0; ++i) 
//    for (int j=ly0; j<=my0; ++j)
//      AzStream << i << " " << j << " " << Az(i,j) << "\n";
//  AzStream.close();
    
//  for (int i=lx0; i<=mx0; ++i) 
//    for (int j=ly0; j<=my0; ++j) 
//      Az(i,j) *= -1;
           
  /* *************************************
   * ... resulting in Bx and By
   */

//  std::cerr << "dx1="<<dx[1]<<"\n";

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      Bx(i,j) = (Az(i,j+1) - Az(i,j-1)) / (2*dx[1]);

//  std::ofstream BxStream("BX.out");
//  for (int i=lx0; i<=mx0; ++i) 
//    for (int j=ly0; j<=my0; ++j)
//      BxStream << i << " " << j << " " << Bx(i,j) << "\n";
//  AzStream.close();

  bound.ScalarFieldReduce(Bx);
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      By(i,j) = (Az(i-1,j) - Az(i+1,j)) / (2*dx[0]);

  bound.ScalarFieldReduce(By);

  /* *************************************
   * The z--component of  -rot(j).
   */

  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      In(i,j) = (jy(i+1,j) - jy(i-1,j)) / (2*dx[0])
		     - (jx(i,j+1) - jx(i,j-1)) / (2*dx[1]);

  In.setParity(ScalarField::EvenParity);
  In.setComponent(ScalarField::ZComponent); 
  bound.ScalarFieldReduce(In);

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

  pois->solve(Bz,In,bound.getNumBoundary(Bz));
    
//  double Bzsum = 0;
//  for (int i=lx1; i<=mx1; ++i) 
//    for (int j=ly1; j<=my1; ++j) {
//      Bzsum += Bz(i,j);
//  }
//  Bzsum /= double( (mx0-lx1)*(my0-ly1) );
//  
//  for (int i=lx0; i<=mx0; ++i) 
//    for (int j=ly0; j<=my0; ++j) {
//      Bz(i,j) -= Bzsum;
//  }
  
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
      Lambda(i,j) = om2(i,j);
    }
  }
           
  // x-component of transverse electric field
  for (int i=lx1; i<=mx1; ++i) {
    for (int j=ly1; j<=my1; ++j) {
      In(i,j) = -(vxx(i+1,j) - vxx(i-1,j)) / (2*dx[0])    
                -(vxy(i,j+1) - vxy(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
                +om2(i,j)*Ex(i,j)                      /// om2*E
                +sy(i,j)*Bz(i,j)-sz(i,j)*By(i,j);      /// q/m rho <v> x B
    }
  }

  In.setParity(ScalarField::OddParity);
  In.setComponent(ScalarField::XComponent); 
  bound.ScalarFieldReduce(In);

  helmh->solve(Out,In,Lambda,bound.getNumBoundary(Etx));
    
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      Etx(i,j) = -Out(i,j);
    
  // y-component of transverse electric field
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) {
      In(i,j) = -(vxy(i+1,j) - vxy(i-1,j)) / (2*dx[0])    
                -(vyy(i,j+1) - vyy(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
                +om2(i,j)*Ey(i,j)                   /// om2*E
                +sz(i,j)*Bx(i,j)-sx(i,j)*Bz(i,j);    /// q/m rho <v> x B
    }

  In.setParity(ScalarField::OddParity);
  In.setComponent(ScalarField::YComponent); 
  bound.ScalarFieldReduce(In);

  helmh->solve(Out,In,Lambda,bound.getNumBoundary(Ety));
    
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      Ety(i,j) = -Out(i,j);
    }
    
  // z-component of transverse electric field
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) {
      In(i,j) = -(vxz(i+1,j) - vxz(i-1,j)) / (2*dx[0])    
                -(vyz(i,j+1) - vyz(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
                +sx(i,j)*By(i,j)-sy(i,j)*Bx(i,j); /// q/m rho <v> x B
    }
    
  In.setParity(ScalarField::OddParity);
  In.setComponent(ScalarField::ZComponent); 
  bound.ScalarFieldReduce(In);

  helmh->solve(Out,In,Lambda,bound.getNumBoundary(Ez));
    
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      Ez(i,j) = -Out(i,j);
//      Ez(i,j) = 0;
    }
           
  clearDiv(Etx, Ety);

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) {
      Ex(i,j) += Etx(i,j);
      Ey(i,j) += Ety(i,j);
    }
        
  double h=1./((mx1-lx0)*(my1-ly0));

  double sumEx=0, sumEy=0, sumEz=0;

  for(int j = ly1; j <= my1; j++) {
    for(int i = lx1; i <= mx1; i++) {
      sumEx += Ex(i,j);
      sumEy += Ey(i,j);
//      sumEz += Ez(i,j);
    }
  }

  sumEx *= h;
  sumEy *= h;
  sumEz *= h;

  for(int j = ly0; j <= my0; j++) {
    for(int i = lx0; i <= mx0; i++) {
      Ex(i,j) -= sumEx;
      Ey(i,j) -= sumEy;
//      Ez(i,j) -= sumEz;
    }
  }

  /* *************************************
   *  This should be it!
   */

//  cerr << "Done DARWIN\n";

  return false;
    
}

void Darwin::clearDiv(ScalarField &Fx, ScalarField &Fy) {
  int lx0 = LBound[0], lx1 = LBound[0]+1;
  int ly0 = LBound[1], ly1 = LBound[1]+1;
  int mx0 = HBound[0], mx1 = HBound[0]-1;
  int my0 = HBound[1], my1 = HBound[1]-1;
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j)  
      DivEt(i,j) = -(Fx(i+1,j) - Fx(i-1,j) ) /(2*dx[0])
	    -(Fy(i,j+1) - Fy(i,j-1) ) /(2*dx[1]);
     
  const Boundary &bound = Process::instance().getBoundary();
  bound.ScalarFieldReduce(DivEt);

  pois->solve(Theta,DivEt,bound.getNumBoundary(Theta));
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) { 
      Fx(i,j) -= (Theta(i+1,j) - Theta(i-1,j) ) /(2*dx[0]);
      Fy(i,j) -= (Theta(i,j+1) - Theta(i,j-1) ) /(2*dx[1]);
    }
  
  bound.ScalarFieldReduce(Fx);
  bound.ScalarFieldReduce(Fy);
}


ScalarField *Darwin::GetByName(const std::string& name) {
  if ("Ex"==name)
    return &Ex; 
  else if ("Ey"==name)
    return &Ey; 
  else if ("Ez"==name)
    return &Ez; 
  else if ("Bx"==name)
    return &Bx; 
  else if ("By"==name)
    return &By; 
  else 
    return &Bz; 
}

