#include <strstream>
#include <fstream>
#include "vlasov.h"
#include "boundary.h"

/* **************************************************************
 *                  SinglePeriodicBoundary                      *
 ****************************************************************/

void SinglePeriodicBoundary::exchangeX(VlasovDist &field) {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    PositionI Xi;
    VelocityI Vi;

    int mx0=UBound[0], mx1=mx0-1, mx2=mx0-2, mx3=mx0-3;
    int lx0=LBound[0], lx1=lx0+1, lx2=lx0+2, lx3=lx0+3;
    
    for (Xi[1] = LBound[1]; Xi[1] <= UBound[1]; ++Xi[1])
      for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
              field(lx0, Xi[1], Vi[0], Vi[1], Vi[2]) =
                    field(mx3, Xi[1], Vi[0], Vi[1], Vi[2]);
              field(lx1, Xi[1], Vi[0], Vi[1], Vi[2]) =
                    field(mx2, Xi[1], Vi[0], Vi[1], Vi[2]);
              field(mx0, Xi[1], Vi[0], Vi[1], Vi[2]) =
                    field(lx3, Xi[1], Vi[0], Vi[1], Vi[2]);
              field(mx1, Xi[1], Vi[0], Vi[1], Vi[2]) =
                    field(lx2, Xi[1], Vi[0], Vi[1], Vi[2]);
          }
}

void SinglePeriodicBoundary::exchangeY(VlasovDist &field) {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    PositionI Xi;
    VelocityI Vi;

    int my0=UBound[1], my1=my0-1, my2=my0-2, my3=my0-3;
    int ly0=LBound[1], ly1=ly0+1, ly2=ly0+2, ly3=ly0+3;
    
    for (Xi[0] = LBound[0]; Xi[0] <= UBound[0]; ++Xi[0])
      for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
              field(Xi[0], ly0, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], my3, Vi[0], Vi[1], Vi[2]);
              field(Xi[0], ly1, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], my2, Vi[0], Vi[1], Vi[2]);
              field(Xi[0], my0, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], ly3, Vi[0], Vi[1], Vi[2]);
              field(Xi[0], my1, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], ly2, Vi[0], Vi[1], Vi[2]);
          }
}

void SinglePeriodicBoundary::ScalarFieldReduce(ScalarField &field) {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    // only one boundary cell
    int mx0=UBound[0], mx1=mx0-1;
    int lx0=LBound[0], lx1=lx0+1;

    int my0=UBound[1], my1=my0-1;
    int ly0=LBound[1], ly1=ly0+1;
    
    for (int i = lx0; i<=mx0; ++i) {
        field(i,ly0) = field(i,my1);
        field(i,my0) = field(i,ly1);
    }
    for (int j = ly0; j<=my0; ++j) {
        field(lx0,j) = field(mx1,j);
        field(mx0,j) = field(lx1,j);
    }       
}

PhasePositionI &SinglePeriodicBoundary::DistLow() {
    return GlDistLow;
}

PhasePositionI &SinglePeriodicBoundary::DistHigh() {
    return GlDistHigh;
}


/* **************************************************************
 *                 MPIPeriodicSplitXBoundary                    *
 ****************************************************************/

MPIPeriodicSplitXBoundary::MPIPeriodicSplitXBoundary(int argc, char **argv) :
    Low(GlDistLow), High(GlDistHigh) {
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&ComSize);
    
    int periodic = true;
    
    MPI_Cart_create(MPI_COMM_WORLD,1,&ComSize,&periodic,true,&comm); 
    MPI_Comm_rank(comm,&ComRank);
    
    MPI_Cart_coords(comm,ComRank,1,&mycoord); 
    
    MPI_Cart_shift(comm,0,1,&leftcoord,&rightcoord); 
            
    double width = (High[0]-2.)/double(ComSize);
    
    if (ComRank>0) 
        Low[0] = int(width*mycoord)-1;
    
    if (ComRank<(ComSize-1))
        High[0] = int(width*(mycoord+1))+2;
        
    exchSize = 2*(High[1]-Low[1]+1)
                *(High[2]-Low[2]+1)
                *(High[3]-Low[3]+1)
                *(High[4]-Low[4]+1);
                
    sendarr = new double[exchSize];
    recvarr = new double[exchSize];
    
    scalarSize = (GlGridHigh[0]-GlGridLow[0]-1)
                *(GlGridHigh[1]-GlGridLow[1]-1);
    scalarsend = new double[scalarSize];
    scalarrecv = new double[scalarSize];
}

MPIPeriodicSplitXBoundary::~MPIPeriodicSplitXBoundary() {
    MPI_Finalize();
    delete[] sendarr;
    delete[] recvarr;
    delete[] scalarrecv;
    delete[] scalarsend;
}

void MPIPeriodicSplitXBoundary::exchangeX(VlasovDist &field) {
    PositionI Xi;
    VelocityI Vi;

    MPI_Status stat; 
    
    int arr_ind = 0;
    for (Xi[0] = Low[0]+2; Xi[0] <= Low[0]+3; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              sendarr[arr_ind++] = 
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);

    MPI_Sendrecv(sendarr, exchSize, MPI_DOUBLE, leftcoord, 0, 
                 recvarr, exchSize, MPI_DOUBLE, rightcoord, 0, 
                 comm, &stat); 

    arr_ind = 0;
    for (Xi[0] = High[0]-1; Xi[0] <= High[0]; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                = recvarr[arr_ind++]; 
                
    arr_ind = 0;
    for (Xi[0] = High[0]-3; Xi[0] <= High[0]-2; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              sendarr[arr_ind++] = 
                field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                
    MPI_Sendrecv(sendarr, exchSize, MPI_DOUBLE, rightcoord, 0, 
                 recvarr, exchSize, MPI_DOUBLE, leftcoord, 0, 
                 comm, &stat); 

    arr_ind = 0;
    for (Xi[0] = Low[0]; Xi[0] <= Low[0]+1; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                = recvarr[arr_ind++]; 
                
}

void MPIPeriodicSplitXBoundary::exchangeY(VlasovDist &field) {
    const int *UBound = field.getHigh();
    const int *LBound = field.getLow();
    
    PositionI Xi;
    VelocityI Vi;

    int my0=UBound[1], my1=my0-1, my2=my0-2, my3=my0-3;
    int ly0=LBound[1], ly1=ly0+1, ly2=ly0+2, ly3=ly0+3;
    
    for (Xi[0] = LBound[0]; Xi[0] <= UBound[0]; ++Xi[0])
      for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
              field(Xi[0], ly0, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], my3, Vi[0], Vi[1], Vi[2]);
              field(Xi[0], ly1, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], my2, Vi[0], Vi[1], Vi[2]);
              field(Xi[0], my0, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], ly3, Vi[0], Vi[1], Vi[2]);
              field(Xi[0], my1, Vi[0], Vi[1], Vi[2]) =
                    field(Xi[0], ly2, Vi[0], Vi[1], Vi[2]);
          }
}

void MPIPeriodicSplitXBoundary::ScalarFieldReduce(ScalarField &field) {
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

    for (int i = lx0; i<=mx0; ++i) {
        field(i,ly0) = field(i,my1);
        field(i,my0) = field(i,ly1);
    }
    for (int j = ly0; j<=my0; ++j) {
        field(lx0,j) = field(mx1,j);
        field(mx0,j) = field(lx1,j);
    }       
}

PhasePositionI &MPIPeriodicSplitXBoundary::DistLow() {
    return Low;
}

PhasePositionI &MPIPeriodicSplitXBoundary::DistHigh() {
    return High;
}

double MPIPeriodicSplitXBoundary::AvgReduce(double val) {
    double result;
    MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
    return result/double(ComSize);
}

/* **************************************************************
 *                 MPIPeriodicSplitXYBoundary                    *
 ****************************************************************/

MPIPeriodicSplitXYBoundary::MPIPeriodicSplitXYBoundary(int argc, char **argv) :
    Low(GlDistLow), High(GlDistHigh) {
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&ComSize);
    
    int periodic[2] = { true, true};
    int dims[2];
    
    dims[0] = int(sqrt(ComSize+0.01));
    dims[1] = int(ComSize/dims[0]);
    
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periodic,true,&comm); 
    MPI_Comm_rank(comm,&ComRank);
    
    MPI_Cart_coords(comm,ComRank,2,mycoord); 
    
    MPI_Cart_shift(comm,0,1,&leftcoord,&rightcoord); 
    MPI_Cart_shift(comm,1,1,&topcoord,&bottomcoord); 
            
    double width[2];
    width[0] = (High[0]-2.)/double(dims[0]);
    width[1] = (High[1]-2.)/double(dims[1]);
    
    if (mycoord[0]>0) 
        Low[0] = int(width[0]*mycoord[0])-1;
    
    if (mycoord[1]>0) 
        Low[1] = int(width[0]*mycoord[1])-1;
    
    if (mycoord[0]<(dims[0]-1))
        High[0] = int(width[1]*(mycoord[0]+1))+2;

    if (mycoord[1]<(dims[1]-1))
        High[1] = int(width[1]*(mycoord[1]+1))+2;
        
    exchSize[0] = 2*(High[1]-Low[1]+1)
                   *(High[2]-Low[2]+1)
                   *(High[3]-Low[3]+1)
                   *(High[4]-Low[4]+1);
    exchSize[1] = 2*(High[0]-Low[0]+1)
                   *(High[2]-Low[2]+1)
                   *(High[3]-Low[3]+1)
                   *(High[4]-Low[4]+1);
                
    sendarrx = new double[exchSize[0]];
    recvarrx = new double[exchSize[0]];
    sendarry = new double[exchSize[1]];
    recvarry = new double[exchSize[1]];
    
    scalarSize = (GlGridHigh[0]-GlGridLow[0]-1)
                *(GlGridHigh[1]-GlGridLow[1]-1);
    scalarsend = new double[scalarSize];
    scalarrecv = new double[scalarSize];
    
    char fname[20];
    
    ostrstream S(fname,20);
    S << "boundary"<<ComRank<<".dat"<<char(0);
    ofstream O(fname);
    O << "Coord: " << mycoord[0] << " " << mycoord[1] << "\n";
    O << "Dims: " << dims[0] << " " << dims[1] << "\n";
    O << "Low: " << Low << "\n";
    O << "High: " << High << "\n";
    O.close();

}

MPIPeriodicSplitXYBoundary::~MPIPeriodicSplitXYBoundary() {
    MPI_Finalize();
    delete[] sendarrx;
    delete[] recvarrx;
    delete[] sendarry;
    delete[] recvarry;
    delete[] scalarrecv;
    delete[] scalarsend;
}

void MPIPeriodicSplitXYBoundary::exchangeX(VlasovDist &field) {
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

    arr_ind = 0;
    for (Xi[0] = High[0]-1; Xi[0] <= High[0]; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                = recvarrx[arr_ind++]; 
                
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

    arr_ind = 0;
    for (Xi[0] = Low[0]; Xi[0] <= Low[0]+1; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                = recvarrx[arr_ind++]; 
                
}

void MPIPeriodicSplitXYBoundary::exchangeY(VlasovDist &field) {
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

    MPI_Sendrecv(sendarry, exchSize[1], MPI_DOUBLE, topcoord, 0, 
                 recvarry, exchSize[1], MPI_DOUBLE, bottomcoord, 0, 
                 comm, &stat); 

    arr_ind = 0;
    for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
      for (Xi[1] = High[1]-1; Xi[1] <= High[1]; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                = recvarry[arr_ind++]; 
                
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

    arr_ind = 0;
    for (Xi[0] = Low[0]; Xi[0] <= High[0]; ++Xi[0])
      for (Xi[1] = Low[1]; Xi[1] <= Low[1]+1; ++Xi[1])
        for (Vi[0] = Low[2]; Vi[0] <= High[2]; ++Vi[0]) 
          for (Vi[1] = Low[3]; Vi[1] <= High[3]; ++Vi[1]) 
            for (Vi[2] = Low[4]; Vi[2] <= High[4]; ++Vi[2])
              field(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                = recvarry[arr_ind++]; 
                
}

void MPIPeriodicSplitXYBoundary::ScalarFieldReduce(ScalarField &field) {
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

    for (int i = lx0; i<=mx0; ++i) {
        field(i,ly0) = field(i,my1);
        field(i,my0) = field(i,ly1);
    }
    for (int j = ly0; j<=my0; ++j) {
        field(lx0,j) = field(mx1,j);
        field(mx0,j) = field(lx1,j);
    }       
}

double MPIPeriodicSplitXYBoundary::AvgReduce(double val) {
    double result;
    MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
    return result/double(ComSize);
}

PhasePositionI &MPIPeriodicSplitXYBoundary::DistLow() {
    return Low;
}

PhasePositionI &MPIPeriodicSplitXYBoundary::DistHigh() {
    return High;
}


