// -*- C++ -*-
// $Id$

#include "magnetostatic.h"
#include "globals.h"
#include "process.h"

#include <fstream>
// ----------------------------------------------------------------------
// Magnetostatic

void Magnetostatic::AddSpecies(MagnetostaticForce* pS) { 	
    species.push_back(pMagnetostaticForce(pS));
}

void Magnetostatic::Init () 
{
  firststep = true;
  n0 = Parameters::instance().bgDensity();
  
  //  Fields range from 0 to N+1, where the inner points are
  //  1 to N
  
  const Boundary& boundary = Process::instance().getBoundary();
  
  LBound = boundary.scalarLow();
  HBound = boundary.scalarHigh();
  
  dx[0] = Parameters::instance().gridSpace_x();
  dx[1] = Parameters::instance().gridSpace_y();
  
  dt = Parameters::instance().dt();
  
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
  
  Ax.resize(LBound.Data(),HBound.Data());
  Ax.setComponent(ScalarField::XComponent);
  Ax.setParity(ScalarField::OddParity);

  Ay.resize(LBound.Data(),HBound.Data());
  Ay.setComponent(ScalarField::YComponent);
  Ay.setParity(ScalarField::OddParity);

  Az.resize(LBound.Data(),HBound.Data());
  Az.setComponent(ScalarField::ZComponent);
  Az.setParity(ScalarField::OddParity);

  oldjx.resize(LBound.Data(),HBound.Data());
  oldjx.setComponent(ScalarField::XComponent);
  oldjx.setParity(ScalarField::OddParity);

  oldjy.resize(LBound.Data(),HBound.Data());
  oldjy.setComponent(ScalarField::YComponent);
  oldjy.setParity(ScalarField::OddParity);

  oldjz.resize(LBound.Data(),HBound.Data());
  oldjz.setComponent(ScalarField::ZComponent);
  oldjz.setParity(ScalarField::OddParity);

  Ex.resize(LBound.Data(),HBound.Data());
  Ex.setComponent(ScalarField::XComponent);
  Ex.setParity(ScalarField::OddParity);
  
  Ey.resize(LBound.Data(),HBound.Data());
  Ey.setComponent(ScalarField::YComponent);
  Ey.setParity(ScalarField::OddParity);
  
  Ez.resize(LBound.Data(),HBound.Data());
  Ez.setComponent(ScalarField::ZComponent);
  Ez.setParity(ScalarField::OddParity);

    
  Theta.resize(LBound.Data(),HBound.Data());
  Theta.setComponent(ScalarField::ScalarComponent);
  Theta.setParity(ScalarField::EvenParity);
  
  DivF.resize(LBound.Data(),HBound.Data());
  DivF.setComponent(ScalarField::ScalarComponent);
  DivF.setParity(ScalarField::EvenParity);

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

  vxx.clear();
  vxy.clear();
  vxz.clear();
  vyy.clear();
  vyz.clear();
  vzz.clear();

  Pot.clear();
  Ax.clear();
  Ay.clear();
  Az.clear();

  Ex.clear();
  Ey.clear();
  Ez.clear();

  Theta.clear();
  DivF.clear();

  Bx.clear();
  By.clear();
  Bz.clear();

  pois = new Poisson();

  helmh = new Helmholtz();
  
  In.resize(LBound.Data(),HBound.Data());
  Lambda.resize(LBound.Data(),HBound.Data());
  Out.resize(LBound.Data(),HBound.Data());

  std::cout << "Done Magnetostatic: Size=( " << LBound << "),(" << HBound<< ")" << endl;

}

bool Magnetostatic::Execute () {
    
  int i;
  double dF, dF2;
  
  int lx0 = LBound[0], lx1 = LBound[0]+1;
  int ly0 = LBound[1], ly1 = LBound[1]+1;
  int mx0 = HBound[0], mx1 = HBound[0]-1;
  int my0 = HBound[1], my1 = HBound[1]-1;
  
  const Boundary &bound = Process::instance().getBoundary();
 
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
    MagnetostaticForce* pS = species[s];

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

  // This is the same for every component
  for (int i=lx0; i<=mx0; ++i) {
    for (int j=ly0; j<=my0; ++j) {
      Lambda(i,j) = om2(i,j);
    }
  }

//  std::cerr << "Ex\n";
  for (int j=ly1; j<=my1; ++j) 
    for (int i=lx1; i<=mx1; ++i) { 
      In(i,j) = -(den(i+1,j) - den(i-1,j)) / (2*dx[0])
              +(vxx(i+1,j) - vxx(i-1,j)) / (2*dx[0])    
                +(vxy(i,j+1) - vxy(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
                -sy(i,j)*Bz(i,j)+sz(i,j)*By(i,j);      /// q/m rho <v> x B
    }
    
  In.setParity(ScalarField::OddParity);
  In.setComponent(ScalarField::XComponent); 
  bound.ScalarFieldReduce(In);

  helmh->solve(Ex,In,Lambda,bound.getNumBoundary(Ex));

//  std::cerr << "Ey\n";
  for (int j=ly1; j<=my1; ++j) 
    for (int i=lx1; i<=mx1; ++i) { 
      In(i,j) = -(den(i,j+1) - den(i,j-1)) / (2*dx[1])
              +(vxy(i+1,j) - vxy(i-1,j)) / (2*dx[0])    
                +(vyy(i,j+1) - vyy(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
                -sz(i,j)*Bx(i,j)+sx(i,j)*Bz(i,j);    /// q/m rho <v> x B
    }
    
  In.setParity(ScalarField::OddParity);
  In.setComponent(ScalarField::YComponent); 
  bound.ScalarFieldReduce(In);

  helmh->solve(Ey,In,Lambda,bound.getNumBoundary(Ey));

//  std::cerr << "Ez\n";
  for (int j=ly1; j<=my1; ++j) 
    for (int i=lx1; i<=mx1; ++i) { 
      In(i,j) = +(vxz(i+1,j) - vxz(i-1,j)) / (2*dx[0])    
                +(vyz(i,j+1) - vyz(i,j-1)) / (2*dx[1]) /// -grad (rho <vv>)
                -sx(i,j)*By(i,j)+sy(i,j)*Bx(i,j); /// q/m rho <v> x B
    }
    
  In.setParity(ScalarField::OddParity);
  In.setComponent(ScalarField::ZComponent); 
  bound.ScalarFieldReduce(In);

  helmh->solve(Ez,In,Lambda,bound.getNumBoundary(Ez));
  
  bound.ScalarFieldReduce(Ex);
  bound.ScalarFieldReduce(Ey);
  bound.ScalarFieldReduce(Ez);
  
  //clearDiv(Ax,Ay);
  
  /* *************************************
   * ... resulting in Bx and By
   */

 
//  std::cerr << "Bx\n";
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      In(i,j) = (jz(i,j+1) - jz(i,j-1)) / (2*dx[1]);

  pois->solve(Bx,In, bound.getNumBoundary(Bx));
  bound.ScalarFieldReduce(Bx);
    
//  std::cerr << "By\n";
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      In(i,j) = (jz(i-1,j) - jz(i+1,j)) / (2*dx[0]);

  pois->solve(By,In, bound.getNumBoundary(Bx));
  bound.ScalarFieldReduce(By);

//  std::cerr << "Bz\n";
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      In(i,j) = (jy(i+1,j) - jy(i-1,j)) / (2*dx[0])
              - (jx(i,j+1) - jx(i,j-1)) / (2*dx[1]);

  pois->solve(Bz,In, bound.getNumBoundary(Bz));
  bound.ScalarFieldReduce(Bz);
  


  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) 
    {
      oldjx(i,j) = jx(i,j); 
      oldjy(i,j) = jy(i,j); 
      oldjz(i,j) = jz(i,j); 
    }

 
  const VelocityD &GlB = Parameters::instance().bField();
  for (int i=lx0; i<=mx0; ++i) {
    for (int j=ly0; j<=my0; ++j) {
      Bx(i,j) += GlB[0];
      By(i,j) += GlB[1];
      Bz(i,j) += GlB[2];
    }
  }


  return false;
    
}

void Magnetostatic::clearDiv(ScalarField &Fx, ScalarField &Fy) {
  int lx0 = LBound[0], lx1 = LBound[0]+1;
  int ly0 = LBound[1], ly1 = LBound[1]+1;
  int mx0 = HBound[0], mx1 = HBound[0]-1;
  int my0 = HBound[1], my1 = HBound[1]-1;
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j)  
      DivF(i,j) = -(Fx(i+1,j) - Fx(i-1,j) ) /(2*dx[0])
	    -(Fy(i,j+1) - Fy(i,j-1) ) /(2*dx[1]);
     
  const Boundary &bound = Process::instance().getBoundary();
  bound.ScalarFieldReduce(DivF);

  pois->solve(Theta,DivF,bound.getNumBoundary(Theta));
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) { 
      Fx(i,j) -= (Theta(i+1,j) - Theta(i-1,j) ) /(2*dx[0]);
      Fy(i,j) -= (Theta(i,j+1) - Theta(i,j-1) ) /(2*dx[1]);
    }
  
  bound.ScalarFieldReduce(Fx);
  bound.ScalarFieldReduce(Fy);
}


ScalarField *Magnetostatic::GetByName(const std::string& name) {
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

