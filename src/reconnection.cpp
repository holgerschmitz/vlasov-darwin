// -*- C++ -*-
// $Id$

#include "reconnection.h"

#include <sstream>

#ifndef SINGLE_PROCESSOR

/* **************************************************************
 *                 SimpleReconnectionBoundary                    *
 ****************************************************************/

SimpleReconnectionBoundary::SimpleReconnectionBoundary() :
    MPIPeriodicSplitXBoundary()
{} 


SimpleReconnectionBoundary::SimpleReconnectionBoundary(int argc, char **argv) :
    MPIPeriodicSplitXBoundary(argc,argv)
{}

SimpleReconnectionBoundary::~SimpleReconnectionBoundary() {}

void SimpleReconnectionBoundary::exchangeY(VlasovDist &field) {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    PositionI Xi;
    VelocityI Vi;

    int my0=UBound[1], my1=my0-1, my2=my0-2, my3=my0-3;
    int ly0=LBound[1], ly1=ly0+1, ly2=ly0+2, ly3=ly0+3;
    
    int sumvx = LBound[2]+UBound[2];
    int sumvz = LBound[4]+UBound[4];
    
    for (Xi[0] = LBound[0]; Xi[0] <= UBound[0]; ++Xi[0])
      for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) {
        int vxi = sumvx-Vi[0];
        for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
          int vzi = sumvx-Vi[0];
          for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) { 
              field(Xi[0], ly0, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], ly3, vxi, Vi[1], vzi);
              field(Xi[0], ly1, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], ly2, vxi, Vi[1], vzi);
              field(Xi[0], my0, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], my3, vxi, Vi[1], vzi);
              field(Xi[0], my1, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], my2, vxi, Vi[1], vzi);
          }
        }
      }
}

void SimpleReconnectionBoundary::ScalarFieldCombine(ScalarField &field) const {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    
    // only one boundary cell
    int mx0=UBound[0], mx1=mx0-1;
    int lx0=LBound[0], lx1=lx0+1;

    int my0=UBound[1], my1=my0-1;
    int ly0=LBound[1], ly1=ly0+1;
    
    int arr_ind=0;

    for (int i = lx1; i <= mx1; ++i)
      for (int j = ly1; j <= my1; ++j) {
        scalarsend[arr_ind++] = field(i,j);
      }
    
    MPI_Allreduce(scalarsend, scalarrecv, scalarSize, MPI_DOUBLE, MPI_SUM, comm);

    arr_ind=0;

    for (int i = lx1; i <= mx1; ++i)
      for (int j = ly1; j <= my1; ++j) {
        field(i,j) = scalarrecv[arr_ind++];
      }

    ScalarFieldReduce(field);
}

void SimpleReconnectionBoundary::ScalarFieldReduce(ScalarField &field) const {
    const NumBoundary &bound = getNumBoundary(field);
    bound.apply(field);
}

const NumBoundary& SimpleReconnectionBoundary::getNumBoundary(ScalarField &field) const {
    int compEven = (field.getComponent()==ScalarField::YComponent)?1:-1;
    int even = compEven*field.getParity();
    if (1==even)
      return evenBound;
    else
      return oddBound;
}

#endif // SINGLE_PROCESSOR

//=========================================================================
//===============   VlasovReconnectionInit  ====================================
//=========================================================================

VlasovReconnectionInit::VlasovReconnectionInit() {}

PARAMETERMAP* VlasovReconnectionInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["vz0"] = WParameter(new ParameterValue<double>(&vz0, 0));
  (*pm)["vz1"] = WParameter(new ParameterValue<double>(&vz1, 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N0, 0));
  (*pm)["Ninf"] = WParameter(new ParameterValue<double>(&Ninf, 1));
  (*pm)["lambda"] = WParameter(new ParameterValue<double>(&lambda, 10));
   return pm;
}

VlasovReconnectionInit::~VlasovReconnectionInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void VlasovReconnectionInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  double Ym = 0.5*(GlLow[1]+GlHigh[1]);
  
  double dy = Parameters::instance().gridSpace_y();
  double lambda_norm = lambda/dy; 
  
  PositionI Xi;
  VelocityI Vi;

  VelocityD VTh(v_th[0],v_th[1],v_th[2]);
  VelocityD UStream(0,0,0);
  double Sum=0;
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
      
      double sc = sech( (Xi[1]-Ym)/lambda_norm );
//      double cs = cosh( (Xi[1]-Ym)/lambda_norm );
      
      double N = Ninf + N0*sc*sc;
      
      double vz_pert = -vz1*cos(2*PIl*Xi[0]/Nx)*sin(PIl*Xi[1]/Ny);
      
      UStream[2] = sc*sc*vz0/N + vz_pert;

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
            Sum = Sum + F;
            dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) = F;
          }
 
    }
    std::cerr << "DENSITY SUM " << Sum << std::endl;
  //    std::cout << "Initialized " << Xi[0] << std::endl;
}

//=========================================================================
//===============   VlasovPeriodicReconnectionInit  ====================================
//=========================================================================

VlasovPeriodicReconnectionInit::VlasovPeriodicReconnectionInit() {}

PARAMETERMAP* VlasovPeriodicReconnectionInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["vz0"] = WParameter(new ParameterValue<double>(&vz0, 0));
  (*pm)["vz1"] = WParameter(new ParameterValue<double>(&vz1, 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N0, 0));
  (*pm)["Ninf"] = WParameter(new ParameterValue<double>(&Ninf, 1));
  (*pm)["lambda"] = WParameter(new ParameterValue<double>(&lambda, 10));
   return pm;
}

VlasovPeriodicReconnectionInit::~VlasovPeriodicReconnectionInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void VlasovPeriodicReconnectionInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  std::cerr << "PERIODIC RECONNECTION\n";
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  double Ysheet1 = 0.75*(GlLow[1]+1)+0.25*(GlHigh[1]-2);
  double Ysheet2 = 0.25*(GlLow[1]+1)+0.75*(GlHigh[1]-2);
  
  int GlMid = (GlLow[1] + GlHigh[1]-1)/2;
  
  double dy = Parameters::instance().gridSpace_y();
  double lambda_norm = lambda/dy; 
  
  PositionI Xi;
  VelocityI Vi;

  VelocityD VTh(v_th[0],v_th[1],v_th[2]);
  VelocityD UStream(0,0,0);

  const Boundary& boundary = Process::instance().getBoundary();
  ofstream USTREAMStream;
  if (boundary.master())
  {
    static int cnt=1;
    ostringstream ustrname;
    ustrname << "USTREAM"<<cnt++<<".out";
    USTREAMStream.open(ustrname.str().c_str());
  }
   
  for (Xi[0] = L[0]; Xi[0] <= H[0]; ++Xi[0])
    for (Xi[1] = L[1]; Xi[1] <= H[1]; ++Xi[1]) {
      
      double sc1 = sech( (Xi[1]-Ysheet1)/lambda_norm );
      double sc2 = sech( (Xi[1]-Ysheet2)/lambda_norm );

//      double vz_pert = vz1*cos(2*PIl*Xi[0]/Nx)*sin(2*PIl*Xi[1]/Ny);
      double vz_pert = 0;
      
      double N_pert = vz1*cos(2*PIl*Xi[0]/Nx);
      double N;
      double N1 = (N0-N_pert)*sc1*sc1;
      double N2 = (N0+N_pert)*sc2*sc2;
      
      if ( (Xi[1]>(GlLow[1]+1)) && (Xi[1]<GlMid) )
      {
        N = Ninf + N1;
        UStream[2] =  vz0*sc1*sc1/N + vz_pert;
      } else if (  (Xi[1]>GlMid) && (Xi[1]<(GlHigh[1]-2)) )
      {
        N = Ninf + N2;
        UStream[2] =  -vz0*sc2*sc2/N + vz_pert;
      } else 
      {
        N = Ninf;
        UStream[2] =  vz_pert;
      }   
                  
      if (boundary.master())
        USTREAMStream << Xi[0] << " " << Xi[1] << " " << UStream[2] << "\n";
      
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
        
    }
    if (boundary.master()) USTREAMStream.close();
    
  //    std::cout << "Initialized " << Xi[0] << std::endl;
}


