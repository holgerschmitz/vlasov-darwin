#include "vlasovinit.h"
#include "potential.h"
#include "darwin.h"

// -------------------------------------------------------------------
// VlasovSpecies
// -------------------------------------------------------------------

template<class ForceField>
VlasovSpecies<ForceField>::VlasovSpecies () 
    :  ForceField(),
       RKState(-2) { 
    dt = 0.125;
    Mass = 1;
    Charge = 1;

    NV[0] = GlGridSize_vx;
    NV[1] = GlGridSize_vy;
    NV[2] = GlGridSize_vz;
    
    VRange[0] = GlGridRange_vx;
    VRange[1] = GlGridRange_vy;
    VRange[2] = GlGridRange_vz;
}

template<class ForceField>
VlasovSpecies<ForceField>::~VlasovSpecies () {}


/** @brief Initializes the VlasovSpecies
 *
 *  Additionally calls the Init of the ForceField class template.
 *
 *  @todo Test the normalisation.
 */
template<class ForceField>
void VlasovSpecies<ForceField>::Init() {


    cerr << "INITIALIZING Vlasov Species" << endl;

    pPot = new (typename ForceField::FieldType);
    pPot->Init();
	pPot->AddSpecies(this);

    PhasePositionI Size;
    PositionI Nx = pPot->GetNx();
    
    for (int i=0; i<2; ++i) Size[i]  =Nx[i];
    for (int i=0; i<3; ++i) Size[2+i]=NV[i];
    
    resize(Size);
    
	gRho.resize(Nx.Data());

	dx[0] = GlGridSpace_x;
    dx[1] = GlGridSpace_y;
    
    for (int i=0; i<2; ++i) BoxRange[i]=(Nx[i]-1)*dx[i];
    
    cerr << "Init Base " << BoxRange << " " << VRange << endl;
       
    cerr << "Init Field " << Charge/Mass << " " << VRange << endl;
    ForceField::Init(Charge/Mass, dx[0]/dt, pPot);
    
    cerr << "Creating initializer " << endl;
    VlasovWaveGenInit<
        ForceField
    > TheInit(this);
    cerr << "Executing initializer " << endl;
    
    TheInit.initialise(Distribution, VRange);
    
    cerr << "Finished initializer " << endl;
}


template<class ForceField>
void VlasovSpecies<ForceField>::resize(const PhasePositionI &size_) {
    cerr << "RISIZING Vlasov Advancer Base" << endl;
    int DSize[5];
    for (int i=0; i<5; ++i) DSize[i] = size_[i]+1;
    
    Distribution.resize(DSize);

    int Pi = 0;
    
    for (   int pi = 0; pi<GridSize.length(); ++pi, ++Pi) 
                GridSize[pi] = size_[Pi];
                
    for (   int vi = 0; vi < VelSize.length(); ++vi, ++Pi) 
                VelSize[vi] = size_[Pi];
                
    VelSizeH = (VelSize-1.0) / 2.0;

    cerr << "RISIZING Finite Volume Advancer" << endl;

    TempDist.resize(DSize);
    T1.resize(DSize);
    T2.resize(DSize);
}


template<class ForceField>
void VlasovSpecies<ForceField>::write(ostream &O) {
    const int *DSize = Distribution.getHigh();
    for (int i=0; i<=DSize[0]; ++i)
      for (int j=0; j<=DSize[1]; ++j)
        for (int k=0; k<=DSize[2]; ++k)
          for (int l=0; l<=DSize[3]; ++l) {
            for (int m=0; m<=DSize[4]; ++m) 
                O << i << " " << j << " " << 
                    k << " " << l << " " << 
                    m << " " << Distribution(i,j,k,l,m) << endl;
            O << endl;
          }
            
}


/** @brief Calls the advance method of the Advancer base class.
 * then calls Task::Execute to execute sub-Tasks
 */
template<class ForceField>
bool VlasovSpecies<ForceField>::Execute () {
    cerr << "Potential\n";
    pPot->Execute(dt);
    
    cerr << "Advancing\n";
    advance(dt);
//    write_Distribution(Distribution,"Advanced.out");
    return false;
}

/** @brief Create density with the density() method supplied by
 *  VlasovAdvancerBase and copy the result into a grid, as requested
 *  by the Opar architecture.
 */
template<class ForceField>
void VlasovSpecies<ForceField>::MakeRho () {

//    cerr << "VlasovSpecies::MakeRho" << endl;

    ScalarField SFRho = density();

    int xbnd = SFRho.getHigh()[0];
    int ybnd = SFRho.getHigh()[1];
	for (int i = 0; i <= xbnd; ++i) 
	  for (int j = 0; j <= ybnd; ++j)  {
		gRho(i,j) = SFRho(i,j);
//        cerr << "den " << i << " " << j << " " << gRho(i,j) << endl;
    }

}

template<class ForceField>
ScalarField VlasovSpecies<ForceField>::density() {
    
//    cerr << "density\n";
    const int *S = Distribution.getHigh();
    PositionI GS(GlGridSize_x,GlGridSize_y);
    ScalarField DEN(GS.Data());
    for (int ix=0; ix<GlGridSize_x; ++ix) {
      for (int iy=0; iy<GlGridSize_y; ++iy) {
        DEN(ix,iy) = 0;
        for (int j=0; j<=S[2]; ++j)
          for (int k=0; k<=S[3]; ++k)
            for (int l=0; l<=S[4]; ++l) 
                DEN(ix,iy) += Distribution(ix,iy,j,k,l);
//        cerr << "den " << ix << " " << iy << " " << DEN(ix,iy) << endl;
        }
    }
    
    return DEN;
}

template<class ForceField>
VelocityD VlasovSpecies<ForceField>::getJ(int i, int j) {
    
    const int *S = Distribution.getHigh();
    
    double d;
    VelocityI vi;
    VelocityD J, V; 
    J[0] = 0.;
    J[1] = 0.;
    J[2] = 0.;
    for (vi[0]=0; vi[0]<=S[2]; ++vi[0])
        for (vi[1]=0; vi[1]<=S[3]; ++vi[1])
            for (vi[2]=0; vi[2]<=S[4]; ++vi[2]) {
                V = velocity(vi);
                d = Distribution(i,j,vi[0],vi[1],vi[2]);
                for (int i=0; i<3; ++i)
                    J[i] += V[i]*d;
            }
    
    return J;
}

template<class ForceField>
FixedArray<double,6> VlasovSpecies<ForceField>::getVVTens(int i, int j) {
    
    const int *S = Distribution.getHigh();
    
    VelocityI vi;
    VelocityD V;
    double d;
    FixedArray<double,6> VV;
    
    for (int l=0; l<6; ++l) VV[l]=0;

    for (vi[0]=0; vi[0]<=S[2]; ++vi[0])
        for (vi[1]=0; vi[1]<=S[3]; ++vi[1])
            for (vi[2]=0; vi[2]<=S[4]; ++vi[2]) {
                V = velocity(vi);
                d = Distribution(i,j,vi[0],vi[1],vi[2]);
                VV[0] += d*V[0]*V[0];
                VV[1] += d*V[0]*V[1];
                VV[2] += d*V[0]*V[2];
                VV[3] += d*V[1]*V[1];
                VV[4] += d*V[1]*V[2];
                VV[5] += d*V[2]*V[2];
            }
    
    return VV;
}



template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceSpace_x(double timestep) {
        
    PositionI Xi;
    VelocityI Vi;
    
    int bx = Distribution.getHigh(0);
    int bxp = bx+1;
    int bxm = bx - 1;
    int bxM = bx - 2;
    NumMatrix<double, 1> Flux(&bxp);
    NumMatrix<double, 1> Dj(&bxp);
    
    for (Vi[0] = 1; Vi[0] < Distribution.getHigh(2); ++Vi[0]) {

      double deltax = -timestep*velocity(Vi)[0]/dx[0];
      int deltaI = int(floor(deltax)+1);        // i + deltaI - 1/2 <= X_(i+1/2) < i + deltaI + 1/2
      double alpha = deltax - deltaI + 1;       // 0 <= alpha < 1

      for (Xi[1] = 1; Xi[1] < Distribution.getHigh(1); ++Xi[1])
        for (Vi[1] = 1; Vi[1] < Distribution.getHigh(3); ++Vi[1]) 
          for (Vi[2] = 1; Vi[2] < Distribution.getHigh(4); ++Vi[2]) {
            
             
            for (Xi[0] = 1; Xi[0] < bx; ++Xi[0]) {
                // go deltaI right and restrict to 1 .. bx-1
                int j =  ((Xi[0]+deltaI+bxM) % bxm) + 1;
                                
                Flux(Xi[0]) = interpolateX(Xi, Vi, j, alpha, Distribution);
                Dj(Xi[0]) = Distribution(j   ,Xi[1], Vi[0], Vi[1], Vi[2]);
            }
            // wrap around
            Flux(0) = Flux(bxm);
            Flux(bx) = Flux(1);
            
            Dj(0) = Dj(bxm);
            Dj(bx) = Dj(1);
            
            for (Xi[0] = 1; Xi[0] <= bx; ++Xi[0]) {
                Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) = Flux(Xi[0]-1) - Flux(Xi[0]) + Dj(Xi[0]);
            }
            Distribution(0 ,Xi[1] , Vi[0], Vi[1], Vi[2]) = Distribution(bxm ,Xi[1] , Vi[0], Vi[1], Vi[2]);
        }
    }
    
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceSpace_y(double timestep) {
        
    PositionI Xi;
    VelocityI Vi;
    
    int bx = Distribution.getHigh(1);
    int bxp = bx+1;
    int bxm = bx - 1;
    int bxM = bx - 2;
    NumMatrix<double, 1> Flux(&bxp);
    NumMatrix<double, 1> Dj(&bxp);
    
//      std::cerr << "Advancing " << Xi[0] << std::endl;
    for (Vi[1] = 1; Vi[1] < Distribution.getHigh(3); ++Vi[1]) {

      double deltax = -timestep*velocity(Vi)[1]/dx[1];
      int deltaI = int(floor(deltax)+1);        // i + deltaI - 1/2 <= X_(i+1/2) < i + deltaI + 1/2
      double alpha = deltax - deltaI + 1;       // 0 <= alpha < 1

//      std::cout << "Delta x = " << deltax << " " << velocity(Vi)[0] 
//                << " " << deltaI << " " << alpha << std::endl; 
      for (Xi[0] = 1; Xi[0] < Distribution.getHigh(0); ++Xi[0])
        for (Vi[0] = 1; Vi[0] < Distribution.getHigh(2); ++Vi[0]) 
          for (Vi[2] = 1; Vi[2] < Distribution.getHigh(4); ++Vi[2]) {
            
             
            for (Xi[1] = 1; Xi[1] < bx; ++Xi[1]) {
                // go deltaI right and restrict to 1 .. bx-1
                int j =  ((Xi[1]+deltaI+bxM) % bxm) + 1;
                
//                std::cout << "  -  " << Xi[0] << " " << j << " " 
//                    << jl << " " << jr << " " << std::endl;
                
                Flux(Xi[1]) = interpolateY(Xi, Vi, j, alpha, Distribution);
                Dj(Xi[1]) = Distribution(Xi[0], j , Vi[0], Vi[1], Vi[2]);
            }
            // wrap around
            Flux(0) = Flux(bxm);
            Flux(bx) = Flux(1);
            
            Dj(0) = Dj(bxm);
            Dj(bx) = Dj(1);
            
            for (Xi[1] = 1; Xi[1] <= bx; ++Xi[1]) {
                Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) = Flux(Xi[1]-1) - Flux(Xi[1]) + Dj(Xi[1]);
            }
            Distribution(Xi[0], 0, Vi[0], Vi[1], Vi[2]) = Distribution(Xi[0], bxm , Vi[0], Vi[1], Vi[2]);
        }
    }
    
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceVel_x(double timestep) {
        
    PositionI Xi;
    VelocityI Vi;
    
    int bvx = Distribution.getHigh(2);
    int bvxp = bvx+1;
    int bvxm = bvx - 1;
    int bvxM = bvx - 1;
    NumMatrix<double, 1> Flux(&bvxp);
    NumMatrix<double, 1> Dj(&bvxp);
    
    int j_old = -1;
    
    for (Xi[0] = 1; Xi[0] < Distribution.getHigh(0); ++Xi[0]) 
      for (Xi[1] = 1; Xi[1] < Distribution.getHigh(1); ++Xi[1]) {
      
        for (Vi[1] = 0; Vi[1] <= Distribution.getHigh(3); ++Vi[1]) 
          for (Vi[2] = 0; Vi[2] <= Distribution.getHigh(4); ++Vi[2]) {
            
             
            for (Vi[0] = 0; Vi[0] < bvx; ++Vi[0]) {

                VelocityD Vel = velocity(Vi);
                VelocityD F = Force(Xi,Vel);
                double deltavx = -timestep*F[0]/deltaVx();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                int j = Vi[0]+deltaI;
                // go deltaI right 
                
                if ((j>=bound_minus()) && (j<=bvx-bound_plus())) {
                    Flux(Vi[0]) = interpolateVx(Xi, Vi, j, alpha, Distribution);
                    if (j_old==-1)
                        Dj(Vi[0]) = Distribution(Xi[0], Xi[1], j, Vi[1], Vi[2]);
                    else {
                        Dj(Vi[0]) = 0;
                        for (int jj=j_old+1; jj<=j; ++jj)
                            Dj(Vi[0]) += Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);
                        j_old=j;
                    }
                } else {
                    Flux(Vi[0]) = 0;
                    Dj(Vi[0]) = 0;
                }
            }
            
            for (Vi[0] = 1; Vi[0] < bvx; ++Vi[0]) {
                Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) 
                    = Flux(Vi[0]-1) - Flux(Vi[0]) + Dj(Vi[0]);
            }

        }
    }
    
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceVel_y(double timestep) {
        
    PositionI Xi;
    VelocityI Vi;
    
    int bvx = Distribution.getHigh(3);
    int bvxp = bvx+1;
    int bvxm = bvx - 1;
    int bvxM = bvx - 1;
    NumMatrix<double, 1> Flux(&bvxp);
    NumMatrix<double, 1> Dj(&bvxp);
    
    int j_old = -1;
    
    for (Xi[0] = 1; Xi[0] < Distribution.getHigh(0); ++Xi[0]) 
      for (Xi[1] = 1; Xi[1] < Distribution.getHigh(1); ++Xi[1]) {
      
        for (Vi[0] = 0; Vi[0] <= Distribution.getHigh(2); ++Vi[0]) 
          for (Vi[2] = 0; Vi[2] <= Distribution.getHigh(4); ++Vi[2]) {
            
             
            for (Vi[1] = 0; Vi[1] < bvx; ++Vi[1]) {
                VelocityD Vel = velocity(Vi);
                VelocityD F = Force(Xi,Vel);
                double deltavx = -timestep*F[1]/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = Vi[1]+deltaI;
                
                if ((j>=bound_minus()) && (j<=bvx-bound_plus())) {
                    Flux(Vi[1]) = interpolateVy(Xi, Vi, j, alpha, Distribution);
                    Dj(Vi[1]) = Distribution(Xi[0], Xi[1], Vi[0], j, Vi[2]);
                    if (j_old==-1)
                        Dj(Vi[1]) = Distribution(Xi[0], Xi[1], Vi[0], j, Vi[2]);
                    else {
                        Dj(Vi[1]) = 0;
                        for (int jj=j_old+1; jj<=j; ++jj)
                            Dj(Vi[1]) = Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);
                        j_old=j;
                    }
                } else {
                    Flux(Vi[1]) = 0;
                    Dj(Vi[1]) = 0;
                }
            }
            
            for (Vi[1] = 1; Vi[1] < bvx; ++Vi[1]) {
                Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) = Flux(Vi[1]-1) - Flux(Vi[1]) + Dj(Vi[1]);
            }

        }
    }
    
}


template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceVel_z(double timestep) {
        
    PositionI Xi;
    VelocityI Vi;
    
    int bvx = Distribution.getHigh(4);
    int bvxp = bvx+1;
    int bvxm = bvx - 1;
    int bvxM = bvx - 1;
    NumMatrix<double, 1> Flux(&bvxp);
    NumMatrix<double, 1> Dj(&bvxp);
    
    int j_old = -1;
    
    for (Xi[0] = 1; Xi[0] < Distribution.getHigh(0); ++Xi[0]) 
      for (Xi[1] = 1; Xi[1] < Distribution.getHigh(1); ++Xi[1]) {
      
        for (Vi[0] = 0; Vi[0] <= Distribution.getHigh(2); ++Vi[0]) 
          for (Vi[1] = 0; Vi[1] <= Distribution.getHigh(3); ++Vi[1]) {
            
             
            for (Vi[2] = 0; Vi[2] < bvx; ++Vi[2]) {
                VelocityD Vel = velocity(Vi);
                VelocityD F = Force(Xi,Vel);
                double deltavx = -timestep*F[2]/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = Vi[2]+deltaI;
                
                if ((j>=bound_minus()) && (j<=bvx-bound_plus())) {
                    Flux(Vi[2]) = interpolateVz(Xi, Vi, j, alpha, Distribution);
                    Dj(Vi[2]) = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], j);
                    if (j_old==-1)
                        Dj(Vi[2]) = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], j);
                    else {
                        Dj(Vi[2]) = 0;
                        for (int jj=j_old+1; jj<=j; ++jj)
                            Dj(Vi[2]) = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], jj);
                        j_old=j;
                    }
                } else {
                    Flux(Vi[2]) = 0;
                    Dj(Vi[2]) = 0;
                }
            }
            
            for (Vi[2] = 1; Vi[2] < bvx; ++Vi[2]) {
                Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) = Flux(Vi[2]-1) - Flux(Vi[2]) + Dj(Vi[2]);
            }

        }
    }
    
}



template<class ForceField>
void VlasovSpecies<ForceField>::advance(double timestep) {
    const int *DSize = Distribution.getHigh();

    switch (RKState) {
        case -2: RKState = -1;
                 break;
        case -1: RKState = 0;
                 advanceStepA(timestep); 
                 break;
        case 0: advanceStepB(T1, timestep);         // T1 = k1
        
                for (int i=0; i<=DSize[0]; ++i)
                  for (int j=0; j<=DSize[1]; ++j)
                    for (int k=0; k<=DSize[2]; ++k)
                      for (int l=0; l<=DSize[3]; ++l)
                        for (int m=0; m<=DSize[4]; ++m)
                          Distribution(i,j,k,l,m) 
                              = Distribution(i,j,k,l,m) + 0.5*T1(i,j,k,l,m); 
                        // D = y_n + 1/2*k1
                            
                advanceStepA(timestep); 
                RKState = 1;
                break;
        case 1: advanceStepB(T2, timestep);         // T2 = k2
        
                for (int i=0; i<=DSize[0]; ++i)
                  for (int j=0; j<=DSize[1]; ++j)
                    for (int k=0; k<=DSize[2]; ++k)
                      for (int l=0; l<=DSize[3]; ++l)
                        for (int m=0; m<=DSize[4]; ++m)
                        Distribution(i,j,k,l,m) 
                            = Distribution(i,j,k,l,m) - 1.5*T1(i,j,k,l,m) + 2*T2(i,j,k,l,m);
                        // D = y_n - k1 + 2*k2
                        
                advanceStepA(timestep);
                
                for (int i=0; i<=DSize[0]; ++i)
                  for (int j=0; j<=DSize[1]; ++j)
                    for (int k=0; k<=DSize[2]; ++k)
                      for (int l=0; l<=DSize[3]; ++l)
                        for (int m=0; m<=DSize[4]; ++m)
                        T1(i,j,k,l,m) = (7./6.)*T1(i,j,k,l,m) - (4./3.)*T2(i,j,k,l,m);
                        
                RKState = 2;
                break;
        case 2: advanceStepB(T2, timestep);      // T2 = 1/6*k3
                for (int i=0; i<=DSize[0]; ++i)
                  for (int j=0; j<=DSize[1]; ++j)
                    for (int k=0; k<=DSize[2]; ++k)
                      for (int l=0; l<=DSize[3]; ++l)
                        for (int m=0; m<=DSize[4]; ++m)
                          Distribution(i,j,k,l,m) 
                              = Distribution(i,j,k,l,m) + T1(i,j,k,l,m) + (1/6.)*T2(i,j,k,l,m);          // D = y_n + 1/6*k1 + 2/3*k2 + 1/6*k3
                advanceStepA(timestep); 
                RKState = 0;
                break;
    }
}



template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepA(double timestep) {
        
    const int *DSize = Distribution.getHigh();
    for (int i=0; i<=DSize[0]; ++i)
      for (int j=0; j<=DSize[1]; ++j)
        for (int k=0; k<=DSize[2]; ++k)
          for (int l=0; l<=DSize[3]; ++l)
            for (int m=0; m<=DSize[4]; ++m) 
              TempDist(i,j,k,l,m) = Distribution(i,j,k,l,m);

    InterpolationInitStep(Distribution);
    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepB(VlasovDist &Out, double timestep) {

    advanceVel_x(timestep);
    advanceVel_y(timestep);
    advanceVel_z(timestep);

    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
    
    const int *DSize = Distribution.getHigh();
    for (int i=0; i<=DSize[0]; ++i)
      for (int j=0; j<=DSize[1]; ++j)
        for (int k=0; k<=DSize[2]; ++k)
          for (int l=0; l<=DSize[3]; ++l) 
            for (int m=0; m<=DSize[4]; ++m) {
              Out(i,j,k,l,m) = Distribution(i,j,k,l,m)-TempDist(i,j,k,l,m);
              Distribution(i,j,k,l,m) = TempDist(i,j,k,l,m);
    }
}


template<class ForceField>
void VlasovSpecies<ForceField>::InterpolationInitStep(const VlasovDist &Dist) {
    f_infty = 1;
}

template<class ForceField>
double VlasovSpecies<ForceField>::epsilonLeft(double fj, double fjp) {
    double fdiff = fjp-fj;
    double fexc = 2*(f_infty-fj);
    if (fexc<0) 
        return 0;
    else if (2*fj<fdiff)
        return 2*fj/fdiff;
    else if (fexc < (-fdiff) )
        return fexc/(-fdiff);
    else 
        return 1.;
}

template<class ForceField>
double VlasovSpecies<ForceField>::epsilonRight(double fj, double fjm) {
    double fdiff = fjm-fj;
    double fexc = 2*(f_infty-fj);
    if (fexc<0) 
        return 0;
    else if (2*fj<fdiff)
        return 2*fj/fdiff;
    else if (fexc < (-fdiff) )
        return fexc/(-fdiff);
    else 
        return 1.;
}


template<class ForceField>
double VlasovSpecies<ForceField>::interpolateX(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha,
                const VlasovDist &Dist) {
    double fjp = Dist(j+1,Xi[1],Vi[0],Vi[1],Vi[2]);
    double fj  = Dist(j  ,Xi[1],Vi[0],Vi[1],Vi[2]);
    double fjm = Dist(j-1,Xi[1],Vi[0],Vi[1],Vi[2]);
    
    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double VlasovSpecies<ForceField>::interpolateY(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha,
                const VlasovDist &Dist) {
    double fjp = Dist(Xi[0], j+1,Vi[0],Vi[1],Vi[2]);
    double fj  = Dist(Xi[0], j  ,Vi[0],Vi[1],Vi[2]);
    double fjm = Dist(Xi[0], j-1,Vi[0],Vi[1],Vi[2]);
    
    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double VlasovSpecies<ForceField>::interpolateVx(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha,
                const VlasovDist &Dist) {
    double fjp = Dist(Xi[0],Xi[1],j+1,Vi[1],Vi[2]);
    double fj  = Dist(Xi[0],Xi[1],j  ,Vi[1],Vi[2]);
    double fjm = Dist(Xi[0],Xi[1],j-1,Vi[1],Vi[2]);

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double VlasovSpecies<ForceField>::interpolateVy(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha,
                const VlasovDist &Dist) {
    double fjp = Dist(Xi[0],Xi[1],Vi[0],j+1,Vi[2]);
    double fj  = Dist(Xi[0],Xi[1],Vi[0],j  ,Vi[2]);
    double fjm = Dist(Xi[0],Xi[1],Vi[0],j-1,Vi[2]);

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double VlasovSpecies<ForceField>::interpolateVz(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha,
                const VlasovDist &Dist) {
    double fjp = Dist(Xi[0],Xi[1],Vi[0],Vi[1],j+1);
    double fj  = Dist(Xi[0],Xi[1],Vi[0],Vi[1],j  );
    double fjm = Dist(Xi[0],Xi[1],Vi[0],Vi[1],j-1);

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}
