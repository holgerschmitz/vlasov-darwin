#include "vlasovinit.h"
#include "potential.h"
#include "darwin.h"

// -------------------------------------------------------------------
// VlasovSpecies
// -------------------------------------------------------------------

template<class ForceField>
VlasovSpecies<ForceField>::VlasovSpecies (Boundary *boundary_) 
    :  ForceField(),
       RKState(-2),
       boundary(boundary_),
       t(0) { 
    dt = Gl_dt;
    Mass = Gl_Mass;
    Charge = Gl_Charge;

    
    VRange[0] = GlGridRange_vx;
    VRange[1] = GlGridRange_vy;
    VRange[2] = GlGridRange_vz;
}

template<class ForceField>
VlasovSpecies<ForceField>::~VlasovSpecies () {
    cerr << "Destructing Vlasov Species\n";
    delete boundary;
}


/** @brief Initializes the VlasovSpecies
 *
 *  Additionally calls the Init of the ForceField class template.
 *
 *  @todo Test the normalisation.
 */
template<class ForceField>
void VlasovSpecies<ForceField>::Init() {


    cerr << "INITIALIZING Vlasov Species" << endl;

    pPot = new (typename ForceField::FieldType)(boundary->master());
    pPot->Init();
	pPot->AddSpecies(this);

    PhasePositionI Size;
    PositionI Lowx = pPot->GetLBound();
    PositionI Highx = pPot->GetHBound();
    
    
    resize(boundary->DistLow(),boundary->DistHigh());
    
	gRho.resize(Lowx.Data(),Highx.Data());
	EKin.resize(Lowx.Data(),Highx.Data());

	Jx.resize(Lowx.Data(),Highx.Data());
	Jy.resize(Lowx.Data(),Highx.Data());
	Jz.resize(Lowx.Data(),Highx.Data());

	Vxx.resize(Lowx.Data(),Highx.Data());
	Vxy.resize(Lowx.Data(),Highx.Data());
	Vxz.resize(Lowx.Data(),Highx.Data());
	Vyy.resize(Lowx.Data(),Highx.Data());
	Vyz.resize(Lowx.Data(),Highx.Data());
	Vzz.resize(Lowx.Data(),Highx.Data());

	dx[0] = GlGridSpace_x;
    dx[1] = GlGridSpace_y;
    
    for (int i=0; i<2; ++i) BoxRange[i]=(Highx[i]-Lowx[i]-3)*dx[i];
    
    cerr << "Init Base " << BoxRange << " " << VRange << endl;
       
       
    cerr << "Init Field " << Charge/Mass << " " << VRange << endl;
    ForceField::Init(Charge/Mass, dx[0]/dt, pPot);
    
    cerr << "Creating initializer " << endl;
//    VlasovMaxwellInit<ForceField> TheInit(this);
//    VlasovTwoMaxwellInit<ForceField> TheInit(this);
    VlasovWaveGenInit<ForceField> TheInit(this);
    cerr << "Executing initializer " << endl;
    
    TheInit.initialise(Distribution, VRange);
    
    cerr << "Boundary Exchange " << endl;
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);

    cerr << "DONE INITIALIZING Vlasov Species" << endl;

}

template<class ForceField>
void VlasovSpecies<ForceField>::resize(PhasePositionI &low, PhasePositionI &high) {
    cerr << "RISIZING Vlasov Advancer Base " << low << " " << high << endl;
    
    Distribution.resize(low.Data(),high.Data());

    int Pi = 0;
    
    for (   int pi = 0; pi<GridSize.length(); ++pi, ++Pi) { 
        GridSize[pi] = high[Pi]-low[Pi]-3;
        GridSizeI[pi] = high[Pi]-low[Pi]-3;     
    }         
    for (   int vi = 0; vi < VelSize.length(); ++vi, ++Pi) 
                VelSize[vi] = high[Pi]-low[Pi]; 
                
    VelSizeH = (VelSize-1.0) / 2.0;

    cerr << "RISIZING Finite Volume Advancer" << endl;
    cerr << "VelSize = " << VelSize << "  VelSizeH = " << VelSizeH << "\n";
    TempDist.resize(low.Data(),high.Data());
    T1.resize(low.Data(),high.Data());
    T2.resize(low.Data(),high.Data());
}


template<class ForceField>
void VlasovSpecies<ForceField>::write(ostream &O) {
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();

    for (int i=L[0]; i<=H[0]; ++i)
      for (int j=L[1]; j<=H[1]; ++j)
        for (int k=L[2]; k<=H[2]; ++k)
          for (int l=L[3]; l<=H[3]; ++l) {
            for (int m=L[4]; m<=H[4]; ++m) 
                O << i << " " << j << " " << 
                    k << " " << l << " " << 
                    m << " " << Distribution(i,j,k,l,m) << endl;
            O << endl;
          }
            
}

template<class ForceField>
void VlasovSpecies<ForceField>::writeYVySlice(int x,int vx,int vz, string fname) {
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();
    if (   (x<L[0]+2) || (x>H[0]-2) 
        || (vx<L[2]) || (vx>H[2])
        || (vz<L[4]) || (vz>H[4]))  return;
    
    ofstream Slice(fname.c_str());
    
    for (int j=L[1]+2; j<=H[1]-2; ++j)
        for (int k=L[3]; k<=H[3]; ++k)
            Slice << j << " " << k << " " 
                << Distribution(x,j,vx,k,vz) << endl;
    
    Slice.close();
}

/** @brief Calls the advance method of the Advancer base class.
 * then calls Task::Execute to execute sub-Tasks
 */
template<class ForceField>
void VlasovSpecies<ForceField>::Execute () {
    if ( (t%100) == 0) {
        double err = densityError();
        correctDensityError(err);
    }
    t++;

//    cerr << "Potential("<<t<<")\n";
    pPot->Execute(dt);
    
//    cerr << "Advancing("<<t<<")\n";
    advance(dt);
//    cerr << "Done Advancing("<<t<<")\n";

//    write_Distribution(Distribution,"Advanced.out");
}


/** @brief Create density with the density() method supplied by
 *  VlasovAdvancerBase and copy the result into a grid, as requested
 *  by the Opar architecture.
 */
template<class ForceField>
void VlasovSpecies<ForceField>::MakeRho () {

//    cerr << "VlasovSpecies::MakeRho" << endl;

    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();
    
    gRho.clear();
    
    // The region of the electromagnetic fields is on
    // grid cell smaller in every direction
    // It's important that the density is only calculated on the
    // inner cells, since ScalarFieldReduce simply adds all the densities of
    // the processes.
    for (int ix=L[0]+2; ix<=H[0]-2; ++ix) {
      for (int iy=L[1]+2; iy<=H[1]-2; ++iy) {
        for (int j=L[2]; j<=H[2]; ++j)
          for (int k=L[3]; k<=H[3]; ++k)
            for (int l=L[4]; l<=H[4]; ++l) 
                gRho(ix,iy) += Distribution(ix,iy,j,k,l);
//        cerr << "den " << ix << " " << iy << " " << gRho(ix,iy) << endl;
        }
    }
    
    boundary->ScalarFieldReduce(gRho);
}

template<class ForceField>
void VlasovSpecies<ForceField>::MakeJs () {
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();
    
    VelocityI vi;
    VelocityD V;
    double d;
    
    Jx.clear();
    Jy.clear();
    Jz.clear();

    Vxx.clear();
    Vxy.clear();
    Vxz.clear();
    Vyy.clear();
    Vyz.clear();
    Vzz.clear();

    for (int i=L[0]+2; i<=H[0]-2; ++i) {
      for (int j=L[1]+2; j<=H[1]-2; ++j) {

        for (vi[0]=L[2]; vi[0]<=H[2]; ++vi[0])
          for (vi[1]=L[3]; vi[1]<=H[3]; ++vi[1])
            for (vi[2]=L[4]; vi[2]<=H[4]; ++vi[2]) {
                V = velocity(vi);
                d = Distribution(i,j,vi[0],vi[1],vi[2]);
                Jx(i,j) += V[0]*d;
                Jy(i,j) += V[1]*d;
                Jz(i,j) += V[2]*d;
                
                Vxx(i,j) += V[0]*V[0]*d;
                Vxy(i,j) += V[0]*V[1]*d;
                Vxz(i,j) += V[0]*V[2]*d;
                Vyy(i,j) += V[1]*V[1]*d;
                Vyz(i,j) += V[1]*V[2]*d;
                Vzz(i,j) += V[2]*V[2]*d;
                
            }
      }
    }
    
    boundary->ScalarFieldReduce(Jx);
    boundary->ScalarFieldReduce(Jy);
    boundary->ScalarFieldReduce(Jz);

    boundary->ScalarFieldReduce(Vxx);
    boundary->ScalarFieldReduce(Vxy);
    boundary->ScalarFieldReduce(Vxz);
    boundary->ScalarFieldReduce(Vyy);
    boundary->ScalarFieldReduce(Vyz);
    boundary->ScalarFieldReduce(Vzz);
}

template<class ForceField>
ScalarField &VlasovSpecies<ForceField>::KineticEnergy () {
    const int *L = Distribution.getLow();
    const int *H = Distribution.getHigh();

    VelocityI vi;
    VelocityD V;
    double d;

    EKin.clear();

    for (int i=L[0]+2; i<=H[0]-2; ++i) 
      for (int j=L[1]+2; j<=H[1]-2; ++j)

        for (vi[0]=L[2]; vi[0]<=H[2]; ++vi[0])
          for (vi[1]=L[3]; vi[1]<=H[3]; ++vi[1])
            for (vi[2]=L[4]; vi[2]<=H[4]; ++vi[2]) {
                V = velocity(vi);
                d = Distribution(i,j,vi[0],vi[1],vi[2]);
                EKin(i,j) += (V[0]*V[0]+V[1]*V[1]+V[2]*V[2])*d;
    }
    
    boundary->ScalarFieldReduce(EKin);
    return EKin;
}

template<class ForceField>
double VlasovSpecies<ForceField>::TotalEnergy () {
    ScalarField &kin = KineticEnergy();
    ScalarField &fld = FieldEnergy();
    const int *L = kin.getLow();
    const int *H = kin.getHigh();
    
    double ETot = 0;
    double dxdy = dx[0]*dx[1];
    
    
    for (int i=L[0]+1; i<=H[0]-1; ++i) 
        for (int j=L[1]+1; j<=H[1]-1; ++j) {
            ETot += kin(i,j)+ fld(i,j);
    }
    return ETot;
}

template<class ForceField>
VelocityD VlasovSpecies<ForceField>::getJ(int i, int j) {
    return VelocityD(Jx(i,j), Jy(i,j), Jz(i,j));
}

template<class ForceField>
FixedArray<double,6> VlasovSpecies<ForceField>::getVVTens(int i, int j) {    
    FixedArray<double,6> Result;
    
    Result[0] = Vxx(i,j);
    Result[1] = Vxy(i,j);
    Result[2] = Vxz(i,j);
    Result[3] = Vyy(i,j);
    Result[4] = Vyz(i,j);
    Result[5] = Vzz(i,j);
    
    return Result;
    
}



template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceSpace_x(double timestep) {
        
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lx = LBound[0];
    int bx = UBound[0];
    int bxM = bx - 2;
    
    NumMatrix<double, 1> Flux(&lx,&bx);
    NumMatrix<double, 1> Dj(&lx,&bx);

//    cerr << "Advancing in x\n";
    for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) {
      double deltax = -timestep*velocity(Vi)[0]/dx[0];
      int deltaI = int(floor(deltax)+1);        // i + deltaI - 1/2 <= X_(i+1/2) < i + deltaI + 1/2
                                                // für vernünftige Werte von V ist deltaI = 0 oder 1
      double alpha = deltax - deltaI + 1;       // 0 <= alpha < 1

//      cerr << "Vx: " << Vi[0] << " " << deltax << " " << deltaI << "\n";
      for (Xi[1] = LBound[1]; Xi[1] <= UBound[1]; ++Xi[1])
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            
             
            // Asymmetric since we later use Flux(x) and Flux(x-1)
            for (Xi[0] = LBound[0]+1; Xi[0] <= bxM; ++Xi[0]) {
                // go deltaI right: j is Xi[0] or Xi[0]+1
                int j =  Xi[0]+deltaI;
                                
                if ( (deltaI<0) || (deltaI>1) ) cerr << "X: deltaI out of bounds:" << deltaI << endl;

                Flux(Xi[0]) = interpolateX(Xi, Vi, j, alpha, Distribution);
                Dj(Xi[0]) = Distribution(j   ,Xi[1], Vi[0], Vi[1], Vi[2]);
            }
    
            for (Xi[0] = LBound[0]+2; Xi[0] <= bxM; ++Xi[0]) 
                Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) 
                    = Flux(Xi[0]-1) - Flux(Xi[0]) + Dj(Xi[0]);
          }
    }
    
//    cerr << "Exchanging boundaries... "; 
    boundary->exchangeX(Distribution);
//    cerr << " ...done\n"; 
    
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceSpace_y(double timestep) {
        
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lx = LBound[1];
    int bx = UBound[1];
    int bxM = bx - 2;

    NumMatrix<double, 1> Flux(&lx,&bx);
    NumMatrix<double, 1> Dj(&lx,&bx);
    
//    cerr << "Advancing in y\n";
    for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) {

      double deltax = -timestep*velocity(Vi)[1]/dx[1];
      int deltaI = int(floor(deltax)+1);        // i + deltaI - 1/2 <= X_(i+1/2) < i + deltaI + 1/2
      double alpha = deltax - deltaI + 1;       // 0 <= alpha < 1

      for (Xi[0] = LBound[0]; Xi[0] <= UBound[0]; ++Xi[0])
        for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            
            
            // Asymmetric since we later use Flux(x) and Flux(x-1)
            for (Xi[1] = LBound[1]+1; Xi[1] <= bxM; ++Xi[1]) {
                // go deltaI right: j is Xi[1] or Xi[1]+1
                int j =  Xi[1]+deltaI;
                                
                if ( (deltaI<0) || (deltaI>1) ) cerr << "Y: deltaI out of bounds:" << deltaI << endl;

                Flux(Xi[1]) = interpolateY(Xi, Vi, j, alpha, Distribution);
                Dj(Xi[1]) = Distribution(Xi[0], j , Vi[0], Vi[1], Vi[2]);
            }
            for (Xi[1] = LBound[1]+2; Xi[1] <= bxM; ++Xi[1]) 
                Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) 
                    = Flux(Xi[1]-1) - Flux(Xi[1]) + Dj(Xi[1]);

          }
    }
    
    boundary->exchangeY(Distribution);
    
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceVel_x(double timestep) {
        
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[2];
    int bvx = UBound[2];
    
    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    int j_old; // check initialization inside the loop!!
    
//    cerr << "Advancing in vx\n";
    for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
      for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1]) {
      
        j_old = -1;
        
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            
             
            for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0]) {
//                cerr << " - " << Vi[0] <<  " - " << Vi[1] <<  " - " << Vi[2] <<  " - " 
//                    << Xi[0] <<  " - " << Xi[1] << "\n";
                VelocityD Vel = velocity(Vi);
                VelocityD F = Force(Xi,Vel);
                double deltavx = -timestep*F[0]/deltaVx();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                int j = Vi[0]+deltaI;
                // go deltaI right 

                if ( (deltaI<0) || (deltaI>1) ) cerr << "Vx: deltaI out of bounds:" << deltaI << endl;

//                cerr << "Vx: " << Vi[0] << " " << deltavx << " " << deltaI << "\n";

                Dj(Vi[0]) = 0;
                for (int jj=j_old+1; jj<=j; ++jj)
                    Dj(Vi[0]) += Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);
                j_old=j;
                
                if ((j>=lvx+bound_minus()) && (j<=bvx-bound_plus()))
                    Flux(Vi[0]) = interpolateVx(Xi, Vi, j, alpha, Distribution);
                else
                    Flux(Vi[0]) = 0;

            }
            
            Dj(bvx-1) = 0;
            for (int jj=j_old+1; jj<bvx; ++jj)
                Dj(bvx-1) += Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);

            Flux(lvx) = 0;
            Flux(bvx-1) = 0;

            for (Vi[0] = lvx+1; Vi[0] < bvx; ++Vi[0]) {
                Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) 
                    = Flux(Vi[0]-1) - Flux(Vi[0]) + Dj(Vi[0]);
            }

        }
    }
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);

}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceVel_y(double timestep) {
        
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[3];
    int bvx = UBound[3];

    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    int j_old; // check initialization inside the loop!!
    
//    cerr << "Advancing in vy\n";
    for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
      for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1]) {
       
        j_old = -1;
        
        for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            
             
            for (Vi[1] = lvx; Vi[1] < bvx; ++Vi[1]) {
            
                VelocityD Vel = velocity(Vi);
                VelocityD F = Force(Xi,Vel);
                double deltavx = -timestep*F[1]/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = Vi[1]+deltaI;
                
                if ( (deltaI<0) || (deltaI>1) ) cerr << "Vx: deltaI out of bounds:" << deltaI << endl;

                Dj(Vi[1]) = 0;
                for (int jj=j_old+1; jj<=j; ++jj)
                    Dj(Vi[1]) += Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);
                j_old=j;

                if ((j>=lvx+bound_minus()) && (j<=bvx-bound_plus()))
                    Flux(Vi[1]) = interpolateVy(Xi, Vi, j, alpha, Distribution);
                else
                    Flux(Vi[1]) = 0;

            }
            
            Dj(bvx-1) = 0;
            for (int jj=j_old+1; jj<bvx; ++jj)
                Dj(bvx-1) += Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);

            Flux(lvx) = 0;
            Flux(bvx-1) = 0;

            for (Vi[1] = lvx+1; Vi[1] < bvx; ++Vi[1]) {
                Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) = Flux(Vi[1]-1) - Flux(Vi[1]) + Dj(Vi[1]);
            }

        }
    }
    
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);
}


template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceVel_z(double timestep) {
        
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[4];
    int bvx = UBound[4];

    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    int j_old; // check initialization inside the loop!!
    
//    cerr << "Advancing in vz\n";
    for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
      for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1]) {
      
        j_old = -1;
        
        for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
          for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) {
            
             
            for (Vi[2] = lvx; Vi[2] < bvx; ++Vi[2]) {
            
                VelocityD Vel = velocity(Vi);
                VelocityD F = Force(Xi,Vel);
                double deltavx = -timestep*F[2]/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = Vi[2]+deltaI;
                
                if ( (deltaI<0) || (deltaI>1) ) cerr << "Vx: deltaI out of bounds:" << deltaI << endl;

                Dj(Vi[2]) = 0;
                for (int jj=j_old+1; jj<=j; ++jj)
                    Dj(Vi[2]) += Distribution(Xi[0], Xi[1],Vi[0], Vi[1], jj);
                j_old=j;

                if ((j>=lvx+bound_minus()) && (j<=bvx-bound_plus()))
                    Flux(Vi[2]) = interpolateVz(Xi, Vi, j, alpha, Distribution);
                else
                    Flux(Vi[2]) = 0;

            }
            
            Dj(bvx-1) = 0;
            for (int jj=j_old+1; jj<bvx; ++jj)
                Dj(bvx-1) += Distribution(Xi[0], Xi[1],Vi[0], Vi[1], jj);

            Flux(lvx) = 0;
            Flux(bvx-1) = 0;

            for (Vi[2] = lvx+1; Vi[2] < bvx; ++Vi[2]) {
                Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) = Flux(Vi[2]-1) - Flux(Vi[2]) + Dj(Vi[2]);
            }

        }
    }
    
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);
}



template<class ForceField>
void VlasovSpecies<ForceField>::advance(double timestep) {
    switch (RKState) {
        case -2: RKState = -1;
                 break;
        case -1: RKState = 0;
                advanceStepA(timestep); 
                break;
        case 0: advanceStepFull(timestep); 
                break;
    }
//     const int *UBound = Distribution.getHigh();
//     const int *LBound = Distribution.getLow();
// 
//     switch (RKState) {
//         case -2: RKState = -1;
//                  break;
//         case -1: RKState = 0;
//                 for (int i=LBound[0]; i<=UBound[0]; ++i)
//                   for (int j=LBound[1]; j<=UBound[1]; ++j)
//                     for (int k=LBound[2]; k<=UBound[2]; ++k)
//                       for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                         for (int m=LBound[4]; m<=UBound[4]; ++m)
//                           TempDist(i,j,k,l,m) = Distribution(i,j,k,l,m);
//                 advanceStepA(timestep/3.); 
//                 break;
//         case 0: advanceStepB(timestep/3.);         // Distribution = c0 = d1
// 
//                 for (int i=LBound[0]; i<=UBound[0]; ++i)
//                   for (int j=LBound[1]; j<=UBound[1]; ++j)
//                     for (int k=LBound[2]; k<=UBound[2]; ++k)
//                       for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                         for (int m=LBound[4]; m<=UBound[4]; ++m)
//                           T1(i,j,k,l,m) = Distribution(i,j,k,l,m);
// 
// 
//                 advanceStepA(2*timestep/3.); 
//                 RKState = 1;
//                 break;
//         case 1: advanceStepB(2*timestep/3.);        // Distribution = c1
// 
//                 for (int i=LBound[0]; i<=UBound[0]; ++i)
//                   for (int j=LBound[1]; j<=UBound[1]; ++j)
//                     for (int k=LBound[2]; k<=UBound[2]; ++k)
//                       for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                         for (int m=LBound[4]; m<=UBound[4]; ++m) {
//                           Distribution(i,j,k,l,m) 
//                             = Distribution(i,j,k,l,m) + TempDist(i,j,k,l,m) - T1(i,j,k,l,m);
//                           T2(i,j,k,l,m) = Distribution(i,j,k,l,m); // = d2
//                         }
//                 advanceStepA(3*timestep/4.);
// 
//                 RKState = 2;
//                 break;
//         case 2: advanceStepB(3*timestep/4.);      // Distribution = c2
//                 for (int i=LBound[0]; i<=UBound[0]; ++i)
//                   for (int j=LBound[1]; j<=UBound[1]; ++j)
//                     for (int k=LBound[2]; k<=UBound[2]; ++k)
//                       for (int l=LBound[3]; l<=UBound[3]; ++l) 
//                         for (int m=LBound[4]; m<=UBound[4]; ++m) {
//                           Distribution(i,j,k,l,m) 
//                               =  (1/4.)*TempDist(i,j,k,l,m) 
//                                + (3/4.)*T1(i,j,k,l,m) 
//                                +  Distribution(i,j,k,l,m)
//                                - T2(i,j,k,l,m); 
//                           TempDist(i,j,k,l,m) = Distribution(i,j,k,l,m);
//                         } 
//                 advanceStepA(timestep/3.); 
//                 RKState = 0;
//                 break;
//     }
}



template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepA(double timestep) {
        
    InterpolationInitStep(Distribution);
    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepB(double timestep) {

    advanceVel_x(timestep);
    advanceVel_y(timestep);
    advanceVel_z(timestep);

    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
    
}

template<class ForceField>
void VlasovSpecies<ForceField>
        ::advanceStepFull(double timestep) {

    advanceVel_x(timestep);
    advanceVel_y(timestep);
    advanceVel_z(timestep);

    advanceSpace_x(timestep);
    advanceSpace_y(timestep);
    
}


template<class ForceField>
void VlasovSpecies<ForceField>::InterpolationInitStep(const VlasovDist &Dist) {
    f_infty = 1;
}

template<class ForceField>
double VlasovSpecies<ForceField>::densityError() {
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    double avg = 0;

    for (int i=LBound[0]+2; i<=UBound[0]-2; ++i)
      for (int j=LBound[1]+2; j<=UBound[1]-2; ++j)
        for (int k=LBound[2]; k<=UBound[2]; ++k) 
          for (int l=LBound[3]; l<=UBound[3]; ++l) 
            for (int m=LBound[4]; m<=UBound[4]; ++m) 
              avg += Distribution(i,j,k,l,m);
    
    avg = avg/double((UBound[0]-LBound[0]-3)*(UBound[1]-LBound[1]-3));
    
    cerr << "Partial Error in density: " << avg - 1 << endl;

    avg =  boundary->AvgReduce(avg);
    cerr << "Total Error in density: " << avg - 1 << endl;
    return avg;
}

template<class ForceField>
void VlasovSpecies<ForceField>::correctDensityError(double err) {
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    for (int i=LBound[0]; i<=UBound[0]; ++i)
      for (int j=LBound[1]; j<=UBound[1]; ++j)
        for (int k=LBound[2]; k<=UBound[2]; ++k) 
          for (int l=LBound[3]; l<=UBound[3]; ++l) 
            for (int m=LBound[4]; m<=UBound[4]; ++m) 
              Distribution(i,j,k,l,m) /= err;
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
