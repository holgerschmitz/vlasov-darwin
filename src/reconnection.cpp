// -*- C++ -*-
// $Id$

#include "reconnection.h"
#include "scalarfield.h"

#include <sstream>

#ifndef SINGLE_PROCESSOR

//=========================================================================
//================= SimpleReconnectionBoundary ============================
//=========================================================================

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

    if (mycoord[0]<dims[0]-1)
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
//      std::cerr << "Mirror X - right\n";
      for (Xi[0] = High[0]-1; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(2*High[0]-3-Xi[0], Xi[1], Low[2]+High[2]-Vi[0], Vi[1], Vi[2]); 
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

    if (mycoord[0]>0)
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
//      std::cerr << "Mirror X - left\n";
      for (Xi[0] = Low[0]; Xi[0] <= Low[0]+1; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(2*Low[0]+3-Xi[0], Xi[1], Low[2]+High[2]-Vi[0], Vi[1], Vi[2]); 
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

    if (mycoord[1]<dims[1]-1)
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
//      std::cerr << "Mirror Y - top\n";
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = High[1]-1; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*High[1]-3-Xi[1], Vi[0], Low[3]+High[3]-Vi[1], Vi[2]);
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

    if (mycoord[1]>0)
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
//      std::cerr << "Mirror Y - bottom\n";
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= Low[1]+1; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*Low[1]+3-Xi[1], Vi[0], Low[3]+High[3]-Vi[1], Vi[2]); 
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

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        scalarsend[arr_ind++] = field(i,j);
      }
    
    MPI_Allreduce(scalarsend, scalarrecv, scalarSize, MPI_DOUBLE, MPI_SUM, comm);

    arr_ind=0;

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        field(i,j) = scalarrecv[arr_ind++];
      }

//    ScalarFieldReduce(field);
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
        break;
      case ScalarField::YComponent : 
        if (1==par) return evenYBound; 
        else return evenXBound;
        break;
      case ScalarField::ZComponent : 
        if (1==par) return evenZBound; 
        else return ScalarBound;
        break;
      case ScalarField::ScalarComponent : 
        if (1==par) return ScalarBound; 
        else return evenZBound;
        break;
      case ScalarField::XYComponent : 
        return evenZBound;
        break;
      case ScalarField::XZComponent : 
        return evenYBound;
        break;
      case ScalarField::YZComponent : 
        return evenXBound;
        break;
      default : std::cerr << " !!! Default Boundary !!!\n";
        return ScalarBound;
    }
}

//=========================================================================
//================= GEMReconnectionBoundary ===============================
//=========================================================================

GEMReconnectionBoundary::GEMReconnectionBoundary() :
    MPIPeriodicSplitXYBoundary()
{} 


GEMReconnectionBoundary::GEMReconnectionBoundary(int argc, char **argv) :
    MPIPeriodicSplitXYBoundary(argc,argv)
{}

GEMReconnectionBoundary::~GEMReconnectionBoundary() {}

void GEMReconnectionBoundary::exchangeX(VlasovDist &field) {
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

    if (mycoord[0]<dims[0]-1)
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
//      std::cerr << "Mirror X - right\n";
      for (Xi[0] = High[0]-1; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(2*High[0]-3-Xi[0], Xi[1], 
                    Low[2]+High[2]-Vi[0], Vi[1], Vi[2]); 
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

    if (mycoord[0]>0)
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
//      std::cerr << "Mirror X - left\n";
      for (Xi[0] = Low[0]; Xi[0] <= Low[0]+1; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(2*Low[0]+3-Xi[0], Xi[1], 
                    Low[2]+High[2]-Vi[0], Vi[1], Vi[2]); 
    }
                
}

void GEMReconnectionBoundary::exchangeY(VlasovDist &field) {
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

    if (mycoord[1]<dims[1]-1)
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
//      std::cerr << "Mirror Y - top\n";
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = High[1]-1; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*High[1]-3-Xi[1], 
                    Low[2]+High[2]-Vi[0], Low[3]+High[3]-Vi[1], Low[4]+High[4]-Vi[2]);
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

    if (mycoord[1]>0)
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
//      std::cerr << "Mirror Y - bottom\n";
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= Low[1]+1; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*Low[1]+3-Xi[1], 
                    Vi[0], Low[3]+High[3]-Vi[1], Vi[2]); 
    }
      
}

void GEMReconnectionBoundary::ScalarFieldCombine(ScalarField &field) const {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    
    // only one boundary cell
    int mx0=UBound[0], mx1=mx0-1;
    int lx0=LBound[0], lx1=lx0+1;

    int my0=UBound[1], my1=my0-1;
    int ly0=LBound[1], ly1=ly0+1;
    
    int arr_ind=0;

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        scalarsend[arr_ind++] = field(i,j);
      }
    
    MPI_Allreduce(scalarsend, scalarrecv, scalarSize, MPI_DOUBLE, MPI_SUM, comm);

    arr_ind=0;

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        field(i,j) = scalarrecv[arr_ind++];
      }

//    ScalarFieldReduce(field);
}

void GEMReconnectionBoundary::ScalarFieldReduce(ScalarField &field) const {
    const NumBoundary &bound = getNumBoundary(field);
    bound.apply(field);
}

const NumBoundary& GEMReconnectionBoundary::getNumBoundary(ScalarField &field) const {
    int comp = field.getComponent();
    int par = field.getParity();
    
    switch (comp)
    {
      case ScalarField::XComponent : 
        if (1==par) return evenXBound; 
        else return evenYBound;
        break;
      case ScalarField::YComponent : 
        if (1==par) return evenYBound; 
        else return oddYBound;
        break;
      case ScalarField::ZComponent : 
        if (1==par) return evenZBound; 
        else return oddZBound;
        break;
      case ScalarField::ScalarComponent : 
        if (1==par) return ScalarBound; 
        else return evenZBound;
        break;
      case ScalarField::XYComponent : 
        return evenZBound;
        break;
      case ScalarField::XZComponent : 
        return evenXZBound;
        break;
      case ScalarField::YZComponent : 
        return evenXBound;
        break;
      default : std::cerr << " !!! Default Boundary !!!\n";
        return ScalarBound;
    }
}

//=========================================================================
//================= FullGEMReconnectionBoundary ===========================
//=========================================================================

FullGEMReconnectionBoundary::FullGEMReconnectionBoundary() :
    MPIPeriodicSplitXYBoundary()
{} 


FullGEMReconnectionBoundary::FullGEMReconnectionBoundary(int argc, char **argv) :
    MPIPeriodicSplitXYBoundary(argc,argv)
{}

FullGEMReconnectionBoundary::~FullGEMReconnectionBoundary() {}


void FullGEMReconnectionBoundary::exchangeY(VlasovDist &field) {
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

    if (mycoord[1]<dims[1]-1)
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
//      std::cerr << "Mirror Y - top\n";
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = High[1]-1; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*High[1]-3-Xi[1], 
                    Low[2]+High[2]-Vi[0], Low[3]+High[3]-Vi[1], Low[4]+High[4]-Vi[2]);
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

    if (mycoord[1]>0)
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
//      std::cerr << "Mirror Y - bottom\n";
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= Low[1]+1; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*Low[1]+3-Xi[1], 
                    Low[2]+High[2]-Vi[0], Low[3]+High[3]-Vi[1], Low[4]+High[4]-Vi[2]); 
    }
      
}

void FullGEMReconnectionBoundary::ScalarFieldCombine(ScalarField &field) const {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    
    // only one boundary cell
    int mx0=UBound[0], mx1=mx0-1;
    int lx0=LBound[0], lx1=lx0+1;

    int my0=UBound[1], my1=my0-1;
    int ly0=LBound[1], ly1=ly0+1;
    
    int arr_ind=0;

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        scalarsend[arr_ind++] = field(i,j);
      }
    
    MPI_Allreduce(scalarsend, scalarrecv, scalarSize, MPI_DOUBLE, MPI_SUM, comm);

    arr_ind=0;

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        field(i,j) = scalarrecv[arr_ind++];
      }

//    ScalarFieldReduce(field);
}

void FullGEMReconnectionBoundary::ScalarFieldReduce(ScalarField &field) const {
    const NumBoundary &bound = getNumBoundary(field);
    bound.apply(field);
}

const NumBoundary& FullGEMReconnectionBoundary::getNumBoundary(ScalarField &field) const {
    int comp = field.getComponent();
    int par = field.getParity();
    
    switch (comp)
    {
      case ScalarField::XComponent : 
        if (1==par) return evenXBound; 
        else return evenYBound;
        break;
      case ScalarField::YComponent : 
        if (1==par) return evenYBound; 
        else return oddYBound;
        break;
      case ScalarField::ZComponent : 
        if (1==par) return evenZBound; 
        else return oddZBound;
        break;
      case ScalarField::ScalarComponent : 
        if (1==par) return ScalarBound; 
        else return evenZBound;
        break;
      case ScalarField::XYComponent : 
        return evenZBound;
        break;
      case ScalarField::XZComponent : 
        return evenXZBound;
        break;
      case ScalarField::YZComponent : 
        return evenXBound;
        break;
      default : std::cerr << " !!! Default Boundary !!!\n";
        return ScalarBound;
    }
}

//=========================================================================
//================= OpenGEMReconnectionBoundary ===========================
//=========================================================================

void OpenGEMReconnectionBoundary::EzNumBoundary
  ::apply(NumMatrix<double,2> &u) const
{   
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();

  double Ysheet1 = GlLow[1]+1.5;
  double XSymmHigh = GlHigh[0]-1.5;

  double Ny = GlHigh[1]-GlLow[1]-3;
  double dx = Parameters::instance().gridSpace_x();
  
  int lx = u.getLow(0);
  int mx = u.getHigh(0);
  int ly = u.getLow(1);
  int my = u.getHigh(1);

  BoundNeumann   bl(lx,mx);
  BoundNeumann   br(lx,mx);
  BoundNeumann   bb(ly,my);
  
  BoundDirichlet    bt(ly,my);

  const int fl = bl.factor();
  const int fr = br.factor();
  const int ft = bt.factor();
  const int fb = bb.factor();
  
  
  double inflowAmp = inflowfactor*sin(0.5*PIl*(my-Ysheet1)/Ny);
  
  for(int i = lx; i <= mx; i++) {
    double d = (XSymmHigh-i)*dx/inflowwidth;
    
    u(i,bb.low())  = fb*u(i,bb.lowSrc());
    u(i,bt.high()) = ft*u(i,bt.highSrc()) - 2*inflowAmp*exp(-d*d);
  }

  for(int j = ly; j <= my; j++) {
    u(bl.low(),j)  = fl*u(bl.lowSrc(),j);
    u(br.high(),j) = fr*u(br.highSrc(),j);
  }
}


OpenGEMReconnectionBoundary::OpenGEMReconnectionBoundary() :
    MPIPeriodicSplitXYBoundary()
{} 


OpenGEMReconnectionBoundary::OpenGEMReconnectionBoundary(int argc, char **argv) :
    MPIPeriodicSplitXYBoundary(argc,argv)
{}

OpenGEMReconnectionBoundary::~OpenGEMReconnectionBoundary() {}

PARAMETERMAP* OpenGEMReconnectionBoundary::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["vy_in"] = WParameter(new ParameterValue<double>(&vyinflow, 1));
  (*pm)["h_in"] = WParameter(new ParameterValue<double>(&inflowwidth, 1));
  (*pm)["lobe_bx"] = WParameter(new ParameterValue<double>(&BxLobe, 1));
   return pm;
}

void OpenGEMReconnectionBoundary::exchangeX(VlasovDist &field) {
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

    if (mycoord[0]<dims[0]-1)
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
//      std::cerr << "Mirror X - right\n";
      for (Xi[0] = High[0]-1; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(2*High[0]-3-Xi[0], Xi[1], 
                    Low[2]+High[2]-Vi[0], Vi[1], Vi[2]); 
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

    if (mycoord[0]>0)
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
//      std::cerr << "Zero derivative - left\n";
      for (Xi[0] = Low[0]; Xi[0] <= Low[0]+1; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Low[0]+2, Xi[1], Vi[0], Vi[1], Vi[2]); 
    }
                
}

void OpenGEMReconnectionBoundary::exchangeY(VlasovDist &field) {
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

    if (mycoord[1]<dims[1]-1)
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
//    else
//    {
//      std::cerr << "top remains unchanged\n";
//      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
//        for (Xi[1] = High[1]-1; Xi[1] <= High[1]; ++Xi[1])
//          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
//            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
//              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
//                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
//                  = field(Xi[0], 2*High[1]-3-Xi[1], 
//                    Low[2]+High[2]-Vi[0], Low[3]+High[3]-Vi[1], Low[4]+High[4]-Vi[2]);
//    }
    
                
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

    if (mycoord[1]>0)
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
//      std::cerr << "Mirror Y - bottom\n";
      for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
        for (Xi[1] = Low[1]; Xi[1] <= Low[1]+1; ++Xi[1])
          for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
            for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
              for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                  = field(Xi[0], 2*Low[1]+3-Xi[1], 
                    Vi[0], Low[3]+High[3]-Vi[1], Vi[2]); 
    }
      
}

void OpenGEMReconnectionBoundary::ScalarFieldCombine(ScalarField &field) const {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    
    // only one boundary cell
    int mx0=UBound[0], mx1=mx0-1;
    int lx0=LBound[0], lx1=lx0+1;

    int my0=UBound[1], my1=my0-1;
    int ly0=LBound[1], ly1=ly0+1;
    
    int arr_ind=0;

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        scalarsend[arr_ind++] = field(i,j);
      }
    
    MPI_Allreduce(scalarsend, scalarrecv, scalarSize, MPI_DOUBLE, MPI_SUM, comm);

    arr_ind=0;

    for (int i = lx0; i <= mx0; ++i)
      for (int j = ly0; j <= my0; ++j) {
        field(i,j) = scalarrecv[arr_ind++];
      }

//    ScalarFieldReduce(field);
}

void OpenGEMReconnectionBoundary::ScalarFieldReduce(ScalarField &field) const {
    const NumBoundary &bound = getNumBoundary(field);
    bound.apply(field);
}

const NumBoundary& OpenGEMReconnectionBoundary::getNumBoundary(ScalarField &field) const {
    ScalarField::FieldComponent ftype = field.getFieldType();
    ScalarField::ComponentType comp = field.getComponent();
    
    int par = field.getParity();
    
    switch (ftype)
    {
      case ScalarField::BX :
        bxBound.setOffsetTop(BxLobe);
        return bxBound;
        break;
      case ScalarField::EZ :
        ezBound.setParameters(vyinflow*BxLobe, inflowwidth);
        return ezBound;
        break;
      default:
        switch (comp)
        {
          case ScalarField::XComponent : 
            if (1==par) return evenXBound; 
            else return evenYBound;
            break;
          case ScalarField::YComponent : 
            if (1==par) return evenYBound; 
            else return oddYBound;
            break;
          case ScalarField::ZComponent : 
            if (1==par) return evenZBound; 
            else return oddZBound;
            break;
          case ScalarField::ScalarComponent : 
            if (1==par) return ScalarBound; 
            else return evenZBound;
            break;
          case ScalarField::XYComponent : 
            return evenZBound;
            break;
          case ScalarField::XZComponent : 
            return evenXZBound;
            break;
          case ScalarField::YZComponent : 
            return evenXBound;
            break;
          default : std::cerr << " !!! Default Boundary !!!\n";
            return ScalarBound;
        }
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
      
      double N_pert = vz1*cos(PIl*Xi[0]/Nx);
      double N;
      double N1 = (N0-N_pert)*sc1*sc1;
      double N2 = (N0+N_pert)*sc2*sc2;
//      double N2 = (N0-N_pert)*sc2*sc2;
      
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
//===============   GEMReconnectionInit  ====================================
//=========================================================================

GEMReconnectionInit::GEMReconnectionInit() {}

PARAMETERMAP* GEMReconnectionInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["vz0"] = WParameter(new ParameterValue<double>(&vz0, 0));
  (*pm)["pert"] = WParameter(new ParameterValue<double>(&pert, 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N0, 0));
  (*pm)["Ninf"] = WParameter(new ParameterValue<double>(&Ninf, 1));
  (*pm)["lambda"] = WParameter(new ParameterValue<double>(&lambda, 10));
   return pm;
}

GEMReconnectionInit::~GEMReconnectionInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void GEMReconnectionInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  std::cerr << "GEM RECONNECTION\n";
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  double Ysheet1 = GlLow[1]+1.5;
  double XSymm = GlLow[0]+1.5;
  
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
      
      double N_pert = pert*cos(PIl*(Xi[0]-XSymm)/Nx)
                             *cos(0.5*PIl*(Xi[1]-Ysheet1)/Ny);
      double N = sc1*sc1*exp(N_pert);
      
      UStream[2] =  vz0;
                   
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
//===============   FullGEMReconnectionInit  ==============================
//=========================================================================

FullGEMReconnectionInit::FullGEMReconnectionInit() {}

PARAMETERMAP* FullGEMReconnectionInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["vz0"] = WParameter(new ParameterValue<double>(&vz0, 0));
  (*pm)["pert"] = WParameter(new ParameterValue<double>(&pert, 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N0, 0));
  (*pm)["Ninf"] = WParameter(new ParameterValue<double>(&Ninf, 1));
  (*pm)["lambda"] = WParameter(new ParameterValue<double>(&lambda, 10));
   return pm;
}

FullGEMReconnectionInit::~FullGEMReconnectionInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void FullGEMReconnectionInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  std::cerr << "FULL GEM RECONNECTION\n";
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  double Ysheet1 = (GlLow[1] + GlHigh[1]) / 2.0;
  double XSymm = (GlLow[0] + GlHigh[0]) / 2.0;
    
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
      
      double N_pert = pert*cos(PIl*(Xi[0]-XSymm)/Nx)
                             *cos(0.5*PIl*(Xi[1]-Ysheet1)/Ny);
      double N = sc1*sc1*exp(N_pert);
      
      UStream[2] =  vz0;
                   
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
//===============   OpenGEMReconnectionInit  ==============================
//=========================================================================

OpenGEMReconnectionInit::OpenGEMReconnectionInit() {}

PARAMETERMAP* OpenGEMReconnectionInit::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["vz0"] = WParameter(new ParameterValue<double>(&vz0, 0));
  (*pm)["vy_in"] = WParameter(new ParameterValue<double>(&vyinflow, 0));
  (*pm)["h_in"] = WParameter(new ParameterValue<double>(&inflowwidth, 0));
  (*pm)["pert"] = WParameter(new ParameterValue<double>(&pert, 0));
  (*pm)["Therm_vx"] = WParameter(new ParameterValue<double>(&v_th[0], 1));
  (*pm)["Therm_vy"] = WParameter(new ParameterValue<double>(&v_th[1], 1));
  (*pm)["Therm_vz"] = WParameter(new ParameterValue<double>(&v_th[2], 1));
  (*pm)["N0"] = WParameter(new ParameterValue<double>(&N0, 0));
  (*pm)["Ninf"] = WParameter(new ParameterValue<double>(&Ninf, 1));
  (*pm)["lambda"] = WParameter(new ParameterValue<double>(&lambda, 10));
   return pm;
}

OpenGEMReconnectionInit::~OpenGEMReconnectionInit() {}


/** @brief Do the initialisation
 *   
 *  Iterate through the whole distribution function and assign the appropriate
 *  phase space density to every point in phase space. The Phase space density 
 *  is calculated as a Maxwellian distribution.
 */
void OpenGEMReconnectionInit::initialise(ForceFieldBase *pVlasov) {
  VlasovDist &dist = pVlasov->getDistribution();
  const int *L = dist.getLow();
  const int *H = dist.getHigh();
  
  std::cerr << "OPEN RECONNECTION\n";
  
  PhasePositionI GlLow  = Parameters::instance().distLow();
  PhasePositionI GlHigh = Parameters::instance().distHigh();
  
  double Nx = GlHigh[0]-GlLow[0]-3;
  double Ny = GlHigh[1]-GlLow[1]-3;
  double Ysheet1 = GlLow[1]+1.5;
  double XSymm = GlLow[0]+1.5;
  double XSymmHigh = GlHigh[0]-1.5;
    
  double dx = Parameters::instance().gridSpace_x();
  double dy = Parameters::instance().gridSpace_y();
  double lambda_norm = lambda/dy;
  
  PositionI Xi;
  VelocityI Vi;

  VelocityD VTh(v_th[0],v_th[1],v_th[2]);
  VelocityD UStream(0,0,0);
  VelocityD UStreamBG(0,0,0);

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
      
      double N_pert = pert*cos(PIl*(Xi[0]-XSymm)/Nx)
                             *cos(0.5*PIl*(Xi[1]-Ysheet1)/Ny);
                             
      double N = sc1*sc1*exp(N_pert);
      
      double d = (XSymmHigh-Xi[0])*dx/inflowwidth;
      double inflowAmp = -vyinflow*exp(-d*d)*sin(0.5*PIl*(Xi[1]-Ysheet1)/Ny);
      
      UStreamBG[1] =  inflowAmp;
      UStream[2] =  vz0;
                   
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


            VelocityD vd0m((Vm - UStreamBG)/VTh);

            VelocityD vd0p((Vp - UStreamBG)/VTh);
            
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
//===============   VlasovPeriodicReconnectionInit  =======================
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


