template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceSpace_x(double timestep) {
    
    errmsg = false;
    
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
                                
                if ( (deltaI<0) || (deltaI>1) && !errmsg) {
                   cerr << "X: deltaI out of bounds:" << deltaI << endl;
                   errmsg = true;
                }

                Flux(Xi[0]) = interpolateX(Xi, Vi, j, alpha);
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
void PosFluxCons3rdOrder<ForceField>
        ::advanceSpace_y(double timestep) {
        
    errmsg = false;
    
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
                                
                if ( (deltaI<0) || (deltaI>1) && !errmsg )
                {
                   cerr << "Y: deltaI out of bounds:" << deltaI << endl;
                   errmsg = true;
                }

                Flux(Xi[1]) = interpolateY(Xi, Vi, j, alpha);
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
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_x(double timestep) {
        
    errmsg = false;
    
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
            
             
            for (Vi[0] = lvx; Vi[0] < bvx; ++Vi[0]) {
//                cerr << " - " << Vi[0] <<  " - " << Vi[1] <<  " - " << Vi[2] <<  " - " 
//                    << Xi[0] <<  " - " << Xi[1] << "\n";
                VelocityD Vel = velocity(Vi);
                VelocityD F = Force(Xi,Vel,timestep);
                double deltavx = -F[0]/deltaVx();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                int j = Vi[0]+deltaI;
                // go deltaI right 

                if ( (deltaI<0) || (deltaI>1) && !errmsg ) 
                {
                  cerr << "Vx: deltaI out of bounds:" << deltaI << endl;
                  errmsg = true;
                }

//                cerr << "Vx: " << Vi[0] << " " << deltavx << " " << deltaI << "\n";

                Dj(Vi[0]) = 0;
                for (int jj=j_old+1; jj<=j; ++jj)
                    Dj(Vi[0]) += Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);
                j_old=j;
                
                if ((j>=lvx+bound_minus()) && (j<=bvx-bound_plus()))
                    Flux(Vi[0]) = interpolateVx(Xi, Vi, j, alpha);
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
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_y(double timestep) {
        
    errmsg = false;
    
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
                VelocityD F = Force(Xi,Vel,timestep);
                double deltavx = -F[1]/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = Vi[1]+deltaI;
                
                if ( (deltaI<0) || (deltaI>1) && !errmsg ) 
                {
                  cerr << "Vy: deltaI out of bounds:" << deltaI << endl;
                  errmsg=true;
                }

                Dj(Vi[1]) = 0;
                for (int jj=j_old+1; jj<=j; ++jj)
                    Dj(Vi[1]) += Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);
                j_old=j;

                if ((j>=lvx+bound_minus()) && (j<=bvx-bound_plus()))
                    Flux(Vi[1]) = interpolateVy(Xi, Vi, j, alpha);
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
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_z(double timestep) {
        
    errmsg = false;
    
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
                VelocityD F = Force(Xi,Vel,timestep);
                double deltavx = -F[2]/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = Vi[2]+deltaI;
                
                if ( (deltaI<0) || (deltaI>1) && !errmsg ) 
                {
                  cerr << "Vz: deltaI out of bounds:" << deltaI << endl;
                  errmsg = true;
                }

                Dj(Vi[2]) = 0;
                for (int jj=j_old+1; jj<=j; ++jj)
                    Dj(Vi[2]) += Distribution(Xi[0], Xi[1],Vi[0], Vi[1], jj);
                j_old=j;

                if ((j>=lvx+bound_minus()) && (j<=bvx-bound_plus()))
                    Flux(Vi[2]) = interpolateVz(Xi, Vi, j, alpha);
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
double PosFluxCons3rdOrder<ForceField>::epsilonLeft(double fj, double fjp) {
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
double PosFluxCons3rdOrder<ForceField>::epsilonRight(double fj, double fjm) {
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
double PosFluxCons3rdOrder<ForceField>::interpolateX(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha) {
    double fjp = Distribution(j+1,Xi[1],Vi[0],Vi[1],Vi[2]);
    double fj  = Distribution(j  ,Xi[1],Vi[0],Vi[1],Vi[2]);
    double fjm = Distribution(j-1,Xi[1],Vi[0],Vi[1],Vi[2]);
    
    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double PosFluxCons3rdOrder<ForceField>::interpolateY(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha) {
    double fjp = Distribution(Xi[0], j+1,Vi[0],Vi[1],Vi[2]);
    double fj  = Distribution(Xi[0], j  ,Vi[0],Vi[1],Vi[2]);
    double fjm = Distribution(Xi[0], j-1,Vi[0],Vi[1],Vi[2]);
    
    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double PosFluxCons3rdOrder<ForceField>::interpolateVx(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha) {
    double fjp = Distribution(Xi[0],Xi[1],j+1,Vi[1],Vi[2]);
    double fj  = Distribution(Xi[0],Xi[1],j  ,Vi[1],Vi[2]);
    double fjm = Distribution(Xi[0],Xi[1],j-1,Vi[1],Vi[2]);

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double PosFluxCons3rdOrder<ForceField>::interpolateVy(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha) {
    double fjp = Distribution(Xi[0],Xi[1],Vi[0],j+1,Vi[2]);
    double fj  = Distribution(Xi[0],Xi[1],Vi[0],j  ,Vi[2]);
    double fjm = Distribution(Xi[0],Xi[1],Vi[0],j-1,Vi[2]);

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

template<class ForceField>
double PosFluxCons3rdOrder<ForceField>::interpolateVz(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha) {
    double fjp = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j+1);
    double fj  = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j  );
    double fjm = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j-1);

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector

    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

