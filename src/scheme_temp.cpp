
#include "util.h"

template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceSpace_x(double timestep) {
 //   CheckDensity(Distribution,"XStart");
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
                                
                if ( ((deltaI<0) || (deltaI>1)) && !errmsg) {
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
    
//    CheckDensity(Distribution,"XEnd");
    boundary->exchangeX(Distribution);
//    CheckDensity(Distribution,"XDone");
    
}

template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceSpace_y(double timestep) {
        
//    CheckDensity(Distribution,"YStart");
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
                                
                if ( ((deltaI<0) || (deltaI>1)) && !errmsg )
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
//    CheckDensity(Distribution,"YEnd");
    
    boundary->exchangeY(Distribution);
//    CheckDensity(Distribution,"YDone");
    
}

template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_x(double timestep) {
        
//    CheckDensity(Distribution,"VXStart");
    errmsg = false;
    
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[2];
    int bvx = UBound[2];
    
    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    VelocityD Vflux;
    
    int j_old; 
    double new_f_infty = 0;
    for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
      for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1]) {
             
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            Vflux[1] = Vi[1];
            Vflux[2] = Vi[2];
            j_old = lvx-1;
             
            for (Vi[0] = lvx; Vi[0] < bvx; ++Vi[0]) {
                Vflux[0] = Vi[0]+0.5;
                VelocityD Vel = velocity(Vflux);
                double F = ForceX(Xi,Vel,timestep);
                double deltavx = -F/deltaVx();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                int j = max(min(Vi[0]+deltaI, bvx+1), lvx-1);
                // go deltaI right 

//                if ( ((deltaI<0) || (deltaI>1)) && !errmsg ) 
//                {
//                  cerr << "Vx: deltaI out of bounds:" << deltaI << endl;
//                  errmsg = true;
//                }

                Dj(Vi[0]) = 0;
                for (int jj=j_old+1; jj<=min(j,bvx) ; ++jj)
                  Dj(Vi[0]) += Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);
                j_old=j;
                
                if ( (j>=lvx) && (j<=bvx) )
                {
                  Flux(Vi[0]) = interpolateVx(Xi, Vi, j, alpha);
                }
                else
                  Flux(Vi[0]) = 0;

            }
            
            Dj(bvx) = 0;
            for (int jj=j_old+1; jj<=bvx; ++jj)
              Dj(bvx) += Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);

//            for (Vi[0] = lvx+1; Vi[0] < bvx; ++Vi[0]) {
//                double dn = Flux(Vi[0]-1) - Flux(Vi[0]) + Dj(Vi[0]);
//                if (dn<0) {
//                  
//                  assert( false );
//                  Flux(Vi[0]-1) += dn/2.0;
//                  Flux(Vi[0]) -= dn/2.0;
//                }
////                else if (dn>f_infty) {
////                  Flux(Vi[0]-1) -= (dn-f_infty)/2.0;
////                  Flux(Vi[0]) += (dn-f_infty)/2.0;
////                }
//            }
            
            double &dl = Distribution(Xi[0], Xi[1], lvx, Vi[1], Vi[2]);
            dl  = - Flux(lvx) + Dj(lvx);
    
            if (new_f_infty<dl) new_f_infty=dl;
            
            for (Vi[0] = lvx+1; Vi[0] < bvx; ++Vi[0]) {
                double &dd = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                dd = Flux(Vi[0]-1) - Flux(Vi[0]) + Dj(Vi[0]);
                if (new_f_infty<dd) new_f_infty=dd;
            }
            
            double &dh = Distribution(Xi[0], Xi[1], bvx, Vi[1], Vi[2]);
            dh  = Flux(bvx-1) + Dj(bvx);
            if (new_f_infty<dh) new_f_infty=dh;
              
        }
    }
//    CheckDensity(Distribution,"VXEnd");
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);
//    CheckDensity(Distribution,"VXDone");
    f_infty = boundary->MaxReduce(new_f_infty);
}

template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_y(double timestep) {
        
//    CheckDensity(Distribution,"VYStart");
    errmsg = false;
    
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[3];
    int bvx = UBound[3];

    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    VelocityD Vflux;

    int j_old; 
    double new_f_infty = 0;
    
    for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
      for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1]) {
               
        for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            Vflux[0] = Vi[0];
            Vflux[2] = Vi[2];
            
            j_old = lvx-1;
            
            for (Vi[1] = lvx; Vi[1] < bvx; ++Vi[1]) {
                Vflux[1] = Vi[1]+0.5;
                VelocityD Vel = velocity(Vflux);
                
                double F = ForceY(Xi,Vel,timestep);
                double deltavx = -F/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = max(min(Vi[1]+deltaI, bvx+1), lvx-1);
                
//                if ( ((deltaI<0) || (deltaI>1)) && !errmsg ) 
//                {
//                  cerr << "Vy: deltaI out of bounds:" << deltaI 
//                      << ", " << j << ", " << j_old 
//                      << ", " << lvx << ", " << bvx << endl;
//                  errmsg=true;
//                }

                Dj(Vi[1]) = 0;
                for (int jj=j_old+1; jj<=min(j,bvx); ++jj)
                  Dj(Vi[1]) += Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);
                j_old=j;
                  
                if ( (j>=lvx) && (j<=bvx) )
                {
                  Flux(Vi[1]) = interpolateVy(Xi, Vi, j, alpha);
                  assert( (Flux(Vi[1])>0) || (Flux(Vi[1])<1) );
                }
                else
                  Flux(Vi[1]) = 0;
                

            }
            Dj(bvx) = 0;
            for (int jj=j_old+1; jj<=bvx; ++jj)
              Dj(bvx) += Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);
  
//            for (Vi[1] = lvx+1; Vi[1] < bvx; ++Vi[1]) {
//                double dn = Flux(Vi[1]-1) - Flux(Vi[1]) + Dj(Vi[1]);
//                if (dn<0) {
//                  assert( false );
//                  Flux(Vi[1]-1) += dn/2.0;
//                  Flux(Vi[1]) -= dn/2.0;
//                } 
////                else if (dn>f_infty) {
////                  Flux(Vi[1]-1) -= (dn-f_infty)/2.0;
////                  Flux(Vi[1]) += (dn-f_infty)/2.0;
////                }
//            }

            double &dl = Distribution(Xi[0], Xi[1], Vi[0], lvx, Vi[2]);
            dl  = - Flux(lvx) + Dj(lvx);
            
            if (new_f_infty<dl) new_f_infty=dl;

            for (Vi[1] = lvx+1; Vi[1] < bvx; ++Vi[1]) {
                double &dd = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                dd = Flux(Vi[1]-1) - Flux(Vi[1]) + Dj(Vi[1]);
                if (new_f_infty<dd) new_f_infty=dd;
            }

            double &dh = Distribution(Xi[0], Xi[1], Vi[0], bvx, Vi[2]);
            dh = Flux(bvx-1) + Dj(bvx);
            if (new_f_infty<dh) new_f_infty=dh;
        }
    }
//    CheckDensity(Distribution,"VYEnd");
    
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);
//    CheckDensity(Distribution,"VYDone");
    f_infty = boundary->MaxReduce(new_f_infty);
}


template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_z(double timestep) {
        
//    CheckDensity(Distribution,"VZStart");
    errmsg = false;
    
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[4];
    int bvx = UBound[4];

    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    VelocityD Vflux;

    int j_old; // check initialization inside the loop!!
    double new_f_infty = 0;
    
//    cerr << "Advancing in vz\n";
    for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
      for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1]) {
      
        for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
          for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) {
            Vflux[0] = Vi[0];
            Vflux[1] = Vi[1];
            
            j_old = lvx-1;
             
            for (Vi[2] = lvx; Vi[2] < bvx; ++Vi[2]) {
            
                Vflux[2] = Vi[2]+0.5;
                VelocityD Vel = velocity(Vflux);
                
                double F = ForceZ(Xi,Vel,timestep);
                double deltavx = -F/deltaVy();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = max(min(Vi[2]+deltaI, bvx+1), lvx-1);
                
//                if ( ((deltaI<0) || (deltaI>1)) && !errmsg ) 
//                {
//                  cerr << "Vz: deltaI out of bounds:" << deltaI << endl;
//                  errmsg = true;
//                }

                Dj(Vi[2]) = 0;
                for (int jj=j_old+1; jj<=min(j,bvx) ; ++jj)
                  Dj(Vi[2]) += Distribution(Xi[0], Xi[1],Vi[0], Vi[1], jj);
                j_old=j;
                  
                //Flux(Vi[2]) = interpolateVz(Xi, Vi, j, alpha);
                if ( (j>=lvx) && (j<=bvx) )
                {
                  Flux(Vi[2]) = interpolateVz(Xi, Vi, j, alpha);
                  assert( (Flux(Vi[2])>0) || (Flux(Vi[2])<1) );
                }
                else
                  Flux(Vi[2]) = 0;

            }
            
            Dj(bvx) = 0;
            for (int jj=j_old+1; jj<=bvx; ++jj)
              Dj(bvx) += Distribution(Xi[0], Xi[1], Vi[0], Vi[1], jj);
  
//            for (Vi[2] = lvx+1; Vi[2] < bvx; ++Vi[2]) {
//                double dn = Flux(Vi[2]-1) - Flux(Vi[2]) + Dj(Vi[2]);
//                if (dn<0) {
//                  assert( false );
//                  Flux(Vi[2]-1) += dn/2.0;
//                  Flux(Vi[2]) -= dn/2.0;
//                } 
////                else if (dn>f_infty) {
////                  Flux(Vi[2]-1) -= (dn-f_infty)/2.0;
////                  Flux(Vi[2]) += (dn-f_infty)/2.0;
////                }
//            }

            double &dl = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], lvx);
            dl = - Flux(lvx) + Dj(lvx);
            if (new_f_infty<dl) new_f_infty=dl;
            
            for (Vi[2] = lvx+1; Vi[2] < bvx; ++Vi[2]) {
                double &dd = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                dd = Flux(Vi[2]-1) - Flux(Vi[2]) + Dj(Vi[2]);
                if (new_f_infty<dd) new_f_infty=dd;
            }
            double &dh = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], bvx);
            dh = Flux(bvx-1) + Dj(bvx);
            if (new_f_infty<dh) new_f_infty=dh;
        }
    }
//    CheckDensity(Distribution,"VZEnd");
    
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);
//    CheckDensity(Distribution,"VZDone");

    f_infty = boundary->MaxReduce(new_f_infty);
}
/** @todo The epsilon limiters do not go well with the Runge Kutta integration.
 *  Find a solution that does not involve commenting lines.
 *  @todo Spend some time thinking about and checking the limiter, (velocity dependent force, large timesteps, boundary)
 */
template<class ForceField>
inline double PosFluxCons3rdOrder<ForceField>::epsilonLeft(double fj, double fjp) {
    double fdiff = fjp-fj;
    //    double fexc = 2*(f_infty-fj);
    //  if (fexc<0) 
    //    return 0;

//    if (2*fj<fdiff)
//        return 2*fj/fdiff;
//    double fexc=2*(f_infty-fj);    
//    if (fexc < (-fdiff) )
//        return -fexc/fdiff;
//    else 
//        return 1.;

    double fexc=2*(f_infty-fj);    
    double fdiff_safe = fdiff + double(fdiff==0);

    double result[3];
    result[0] = 1;
    result[1] = 2*fj/fdiff_safe;
    result[2] = -fexc/fdiff_safe;
    
    int res1 = (2*fj<fdiff);
    int res2 = (fexc < (-fdiff) );
    
    return result[res1 + 2*(1-res1)*res2];
}

template<class ForceField>
inline double PosFluxCons3rdOrder<ForceField>::epsilonRight(double fj, double fjm) {
    double fdiff = fjm-fj;
    //    double fexc = 2*(f_infty-fj);
    //    if (fexc<0) 
    //        return 0;
//    if (2*fj<fdiff)
//        return 2*fj/fdiff;
//    double fexc = 2*(f_infty-fj);
//    if (fexc < (-fdiff) )
//        return -fexc/fdiff;
//    else 
//        return 1.;

    double fexc=2*(f_infty-fj);    
    double fdiff_safe = fdiff + double(fdiff==0);

    double result[3];
    result[0] = 1;
    result[1] = 2*fj/fdiff_safe;
    result[2] = -fexc/fdiff_safe;
    
    int res1 = (2*fj<fdiff);
    int res2 = (fexc < (-fdiff) );
    
    return result[res1 + 2*(1-res1)*res2];

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
    double fj  = Distribution(Xi[0],Xi[1],j  ,Vi[1],Vi[2]);
    double fjp, fjm;
    if (j<=Distribution.getLow()[2])
    {
      fjp = Distribution(Xi[0],Xi[1],j+1,Vi[1],Vi[2]);
      fjm = max(0.0,2*fj-fjp);
    } else
    if (j>=Distribution.getHigh()[2])
    {
      fjm = Distribution(Xi[0],Xi[1],j-1,Vi[1],Vi[2]);
      fjp = max(0.0,2*fj-fjm);
    } else
    {
      fjp = Distribution(Xi[0],Xi[1],j+1,Vi[1],Vi[2]);
      fjm = Distribution(Xi[0],Xi[1],j-1,Vi[1],Vi[2]);
    }

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
    double fj  = Distribution(Xi[0],Xi[1],Vi[0],j  ,Vi[2]);
    double fjp, fjm;
    if (j<=Distribution.getLow()[3])
    {
      fjp = Distribution(Xi[0],Xi[1],Vi[0],j+1,Vi[2]);
      fjm = max(0.0,2*fj-fjp);
    } else
    if (j>=Distribution.getHigh()[3])
    {
      fjm = Distribution(Xi[0],Xi[1],Vi[0],j-1,Vi[2]);
      fjp = max(0.0,2*fj-fjm);
    } else
    {
      fjp = Distribution(Xi[0],Xi[1],Vi[0],j+1,Vi[2]);
      fjm = Distribution(Xi[0],Xi[1],Vi[0],j-1,Vi[2]);
    }

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
    double fj  = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j  );
    double fjp, fjm;
    if (j<=Distribution.getLow()[4])
    {
      fjp = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j+1);
      fjm = max(0.0,2*fj-fjp);
    } else
    if (j>=Distribution.getHigh()[4])
    {
      fjm = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j-1);
      fjp = max(0.0,2*fj-fjm);
    } else
    {
      fjp = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j+1);
      fjm = Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j-1);
    }

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector
    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

