// -*- C++ -*-
// $Id$

#include "reconnection.h"
#include "scalarfield.h"

#include <sstream>

#ifndef SINGLE_PROCESSOR

/* **************************************************************
 *                 SimpleReconnectionBoundary                    *
 ****************************************************************/

SimpleReconnectionBoundary::SimpleReconnectionBoundary() :
    MPIPeriodicSplitXYBoundary()
{} 


SimpleReconnectionBoundary::SimpleReconnectionBoundary(int argc, char **argv) :
    MPIPeriodicSplitXYBoundary(argc,argv)
{}

SimpleReconnectionBoundary::~SimpleReconnectionBoundary() {}

void SimpleReconnectionBoundary::exchangeX(VlasovDist &field) {
    PositionI Xi;
    VelocityI Vi;

    MPI_Status stat; 
    
    int arr_ind = 0;
    for (Xi[0] = Low[0]+2; Xi[0] <= Low[0]+3; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              sendarrx[arr_ind++] = 
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);

    MPI_Sendrecv(sendarrx, exchSize[0], MPI_DOUBLE, leftcoord, 0, 
                 recvarrx, exchSize[0], MPI_DOUBLE, rightcoord, 0, 
                 comm, &stat); 

    if (rightcoord != MPI_PROC_NULL)
    {
      arr_ind = 0;
      for (Xi[0] = High[0]-1; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = recvarrx[arr_ind++]; 
    }
    else
    {
      for (Xi[0] = High[0]-1; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(2*High[0]-3-Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]); 
    }
    
                
    arr_ind = 0;
    for (Xi[0] = High[0]-3; Xi[0] <= High[0]-2; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              sendarrx[arr_ind++] = 
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                
    MPI_Sendrecv(sendarrx, exchSize[0], MPI_DOUBLE, rightcoord, 0, 
                 recvarrx, exchSize[0], MPI_DOUBLE, leftcoord, 0, 
                 comm, &stat); 

    if (leftcoord != MPI_PROC_NULL)
    {
      arr_ind = 0;
      for (Xi[0] = Low[0]; Xi[0] <= Low[0]+1; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = recvarrx[arr_ind++]; 
    }
    else
    {
      arr_ind = 0;
      for (Xi[0] = Low[0]; Xi[0] <= Low[0]+1; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(2*Low[0]+3-Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]); 
    }
                
}

void SimpleReconnectionBoundary::exchangeY(VlasovDist &field) {
    PositionI Xi;
    VelocityI Vi;

    MPI_Status stat; 
    
    int arr_ind = 0;
    for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
      for (Xi[1] = Low[1]+2; Xi[1] <= Low[1]+3; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              sendarry[arr_ind++] = 
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);

    MPI_Sendrecv(sendarry, exchSize[1], MPI_DOUBLE, bottomcoord, 0, 
                 recvarry, exchSize[1], MPI_DOUBLE, topcoord, 0, 
                 comm, &stat); 

    if (topcoord != MPI_PROC_NULL)
    {
      arr_ind = 0;
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = High[1]-1; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = recvarry[arr_ind++];
    }
    else
    {
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = High[1]-1; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*High[1]-3-Xi[1], Vi[0], Vi[1], Vi[2]);
    }
    
                
    arr_ind = 0;
    for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
      for (Xi[1] = High[1]-3; Xi[1] <= High[1]-2; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              sendarry[arr_ind++] = 
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                
    MPI_Sendrecv(sendarry, exchSize[1], MPI_DOUBLE, topcoord, 0, 
                 recvarry, exchSize[1], MPI_DOUBLE, bottomcoord, 0, 
                 comm, &stat); 

    if (bottomcoord != MPI_PROC_NULL)
    {
      arr_ind = 0;
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= Low[1]+1; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = recvarry[arr_ind++]; 
    }
    else
    {
      arr_ind = 0;
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= Low[1]+1; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*Low[1]+3-Xi[1], Vi[0], Vi[1], Vi[2]); 
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
    int comp = field.getComponent();
    int par = field.getParity();
    
    switch (comp)
    {
      case ScalarField::XComponent : 
        if (1==par) return evenXBound; 
        else return evenYBound;
      case ScalarField::YComponent : 
        if (1==par) return evenYBound; 
        else evenXBound;
      case ScalarField::ZComponent : 
        if (1==par) return evenZBound; 
        else ScalarBound;
      case ScalarField::ScalarComponent : 
        if (1==par) return ScalarBound; 
        else evenZBound;
      default : return ScalarBound;
    }
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
  
  std::cerr << "SYMMETRIC RECONNECTION\n";
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  double Ysheet1 = GlLow[1]+1.5;
  double Ysheet2 = GlHigh[1]-1.5;
  
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
      
      double N_pert = vz1*cos(2*PIl*Xi[0]/Nx);
      double N;
      double N1 = (N0-N_pert)*sc1*sc1;
      double N2 = (N0+N_pert)*sc2*sc2;
      
      if ( Xi[1]<=GlMid )
      {
        N = N1;
        UStream[2] =  vz0;
          
      } else 
      {
        N = N2;
        UStream[2] =  -vz0;
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


            VelocityD vd0m(Vm/VTh);

            VelocityD vd0p(Vp/VTh);
            
            VelocityD F0;
            
            for (int j=0; j<3; ++j) {
              if (Vi[j]==L[j+2]) {
                F0[j] = 0.5*(erf(vd0p[j]) + 1);
              } else if (Vi[j]==H[j+2]) {
                F0[j] = 0.5*(1 - erf(vd0m[j]));
              } else {
                F0[j] = 0.5*(erf(vd0p[j]) - erf(vd0m[j]));
              }
            }
                
            F = Ninf*F0.product();
            dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) += F;

          }
        
    }
    if (boundary.master()) USTREAMStream.close();
    
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
      
      double N_pert = vz1*cos(2*PIl*Xi[0]/Nx);
      double N;
      double N1 = (N0-N_pert)*sc1*sc1;
      double N2 = (N0+N_pert)*sc2*sc2;
      
      if ( (Xi[1]>(GlLow[1]+1)) && (Xi[1]<GlMid) )
      {
        N = N1;
        UStream[2] =  vz0;
          
      } else if (  (Xi[1]>GlMid) && (Xi[1]<(GlHigh[1]-2)) )
      {
        N = N2;
        UStream[2] =  -vz0;
      } else 
      {
        N = 0;
        UStream[2] =  0;
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


            VelocityD vd0m(Vm/VTh);

            VelocityD vd0p(Vp/VTh);
            
            VelocityD F0;
            
            for (int j=0; j<3; ++j) {
              if (Vi[j]==L[j+2]) {
                F0[j] = 0.5*(erf(vd0p[j]) + 1);
              } else if (Vi[j]==H[j+2]) {
                F0[j] = 0.5*(1 - erf(vd0m[j]));
              } else {
                F0[j] = 0.5*(erf(vd0p[j]) - erf(vd0m[j]));
              }
            }
                
            F = Ninf*F0.product();
            dist(Xi[0],Xi[1],Vi[0],Vi[1],Vi[2]) += F;

          }
        
    }
    if (boundary.master()) USTREAMStream.close();
    
  //    std::cout << "Initialized " << Xi[0] << std::endl;
}


