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
    
  jx.resize(LBound.Data(),HBound.Data());
  jx.setComponent(ScalarField::XComponent);
  jx.setParity(ScalarField::OddParity);
  jy.resize(LBound.Data(),HBound.Data());
  jy.setComponent(ScalarField::YComponent);
  jy.setParity(ScalarField::OddParity);
  jz.resize(LBound.Data(),HBound.Data());
  jz.setComponent(ScalarField::ZComponent);
  jz.setParity(ScalarField::OddParity);
  
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
    
  jx.clear();
  jy.clear();
  jz.clear();

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
  
  In.resize(LBound.Data(),HBound.Data());
  Lambda.resize(LBound.Data(),HBound.Data());
  Out.resize(LBound.Data(),HBound.Data());

  std::cout << "Done Magnetostatic: Size=( " << LBound << "),(" << HBound<< ")" << endl;
}

bool Magnetostatic::Execute () {
    
  int i;
  double dF;
  
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
    
  jx.clear();
  jy.clear();
  jz.clear();
	
  VelocityD jt;
    
  /* *************************************
   *  First ask the species to create the density and the current density
   *  and add them all up.
   */
    
  // iterate through the species
//  std::cerr << "dV="<<dV<<"\n";
  for (int s = species.size() - 1; s >= 0; s--) {
    MagnetostaticForce* pS = species[s];

    dF = pS->getCharge();

    // cerr << "Creating Density\n";
    DistMomentRho *distRho = pS->getDerivedRho();
    DistMomentVelocities *distVel = pS->getDerivedVelocities();
    
    // ... and get it
    ScalarField &rho = distRho->getRho();
    for (int j=ly0; j<=my0; ++j) 
      for (int i=lx0; i<=mx0; ++i) {
	    den(i,j) += dF*rho(i,j);

	    jt = distVel->getJ(i,j);
	    jx(i,j) += dF*jt[0];
	    jy(i,j) += dF*jt[1];
	    jz(i,j) += dF*jt[2];
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
      In(i,j) = (den(i,j)+n0);
    }
    
  pois->solve(Pot,In, bound.getNumBoundary(Pot));

  /* *************************************
   *  Now we want to calculate the magnetic field B.
   *  In two dimensions we calculate the z-component of the
   *  vector potential A_z which leads to the magnetic field
   *  components B_x and B_y.
   *  The component B_z can be calculated directly from rot j_z.
   */

  /* *************************************
   * The x, y and z--component of the vector potential.
   */
  
  pois->solve(Ax,jx, bound.getNumBoundary(Ax));
  pois->solve(Ay,jy, bound.getNumBoundary(Ay));
  pois->solve(Az,jz, bound.getNumBoundary(Az));
  
  clearDiv(Ax,Ay);
  
  /* *************************************
   * ... resulting in Bx and By
   */

//  std::cerr << "dx1="<<dx[1]<<"\n";

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      Bx(i,j) = (Az(i,j+1) - Az(i,j-1)) / (2*dx[1]);

  bound.ScalarFieldReduce(Bx);
    
  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      By(i,j) = (Az(i-1,j) - Az(i+1,j)) / (2*dx[0]);

  bound.ScalarFieldReduce(By);

  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      Bz(i,j) = (Ay(i+1,j) - Ay(i-1,j)) / (2*dx[0])
              - (Ax(i,j+1) - Ax(i,j-1)) / (2*dx[1]);

  bound.ScalarFieldReduce(Bz);
  
  if (firststep)
  {
    firststep=false;
    for (int i=lx0; i<=mx0; ++i) 
      for (int j=ly0; j<=my0; ++j) 
      {
        oldAx(i,j) = Ax(i,j); 
        oldAy(i,j) = Ay(i,j); 
        oldAz(i,j) = Az(i,j); 
      }
  }

  for (int i=lx1; i<=mx1; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      Ex(i,j) = (Pot(i-1,j) - Pot(i+1,j)) / (2*dx[0])
              - (oldAx(i,j) - Ax(i,j)) / dt;

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly1; j<=my1; ++j) 
      Ey(i,j) = (Pot(i,j-1) - Pot(i,j+1)) / (2*dx[1])
              - (oldAy(i,j) - Ay(i,j)) / dt;
              
  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) 
      Ez(i,j) = - (oldAz(i,j) - Az(i,j)) / dt;

  bound.ScalarFieldReduce(Ex);
  bound.ScalarFieldReduce(Ey);
  bound.ScalarFieldReduce(Ez);

  for (int i=lx0; i<=mx0; ++i) 
    for (int j=ly0; j<=my0; ++j) 
    {
      oldAx(i,j) = Ax(i,j); 
      oldAy(i,j) = Ay(i,j); 
      oldAz(i,j) = Az(i,j); 
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

