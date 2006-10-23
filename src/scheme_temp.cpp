#include "util.h" // for function floor(const double&)
//-----------------------------------------------------------------------------
template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceSpace_x(double timestep) {

    errmsg = false;
    //get upper/lower Boundary
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;

    int lx = LBound[0];
    int bx = UBound[0];
    int bxM = bx - 2;

    //create new matrices, with the lower and upper indices given by lx/bx
    NumMatrix<double, 1> Flux(&lx,&bx);
    NumMatrix<double, 1> Dj(&lx,&bx);


    for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) {
	//calculate space step delta x 
      double deltax = -timestep*this->velocity(Vi)[0]/this->dx[0];
	//calculate space step in cells
	// i + deltaI - 1/2 <= X_(i+1/2) < i + deltaI + 1/2
	//for all legitmate values of V: deltaI =  1 or 0
	// 0 <= alpha < 1
      int deltaI = int(floor(deltax)+1);
	//non integer part of deltax
      double alpha = deltax - deltaI + 1;

      for (Xi[1] = LBound[1]+2; Xi[1] <= UBound[1]-2; ++Xi[1])
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            // Asymmetric since we later use Flux(x) and Flux(x-1)
            for (Xi[0] = LBound[0]+1; Xi[0] <= bxM; ++Xi[0]) {
                // go deltaI right: j is Xi[0] or Xi[0]+1
                int j =  Xi[0]+deltaI; 
                if ( ((deltaI<0) || (deltaI>1)) && !errmsg) {
                   cerr << "X: deltaI out of bounds: " << deltaI << endl;
                   errmsg = true;
                }

                Flux(Xi[0]) = interpolateX(Xi, Vi, j, alpha);
                Dj(Xi[0]) = this->Distribution(j,Xi[1], Vi[0], Vi[1], Vi[2]);
            }
            for (Xi[0] = LBound[0]+2; Xi[0] <= bxM; ++Xi[0]) 
                this->Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) 
                    = Flux(Xi[0]-1) - Flux(Xi[0]) + Dj(Xi[0]);
          }
    }
    //exchange new boundaries, in both space directions
    this->boundary->exchangeX(this->Distribution);
    this->boundary->exchangeY(this->Distribution); 
}
//-----------------------------------------------------------------------------
template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceSpace_y(double timestep) {
        
//    CheckDensity(Distribution,"YStart");
    errmsg = false;
    
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lx = LBound[1];
    int bx = UBound[1];
    int bxM = bx - 2;

    NumMatrix<double, 1> Flux(&lx,&bx);
    NumMatrix<double, 1> Dj(&lx,&bx);
    
//    cerr << "Advancing in y\n";
    for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) {

      double deltax = -timestep*this->velocity(Vi)[1]/this->dx[1];
      int deltaI = int(floor(deltax)+1);        // i + deltaI - 1/2 <= X_(i+1/2) < i + deltaI + 1/2
      double alpha = deltax - deltaI + 1;       // 0 <= alpha < 1

      for (Xi[0] = LBound[0]+2; Xi[0] <= UBound[0]-2; ++Xi[0])
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
                Dj(Xi[1]) = this->Distribution(Xi[0], j , Vi[0], Vi[1], Vi[2]);
            }
            for (Xi[1] = LBound[1]+2; Xi[1] <= bxM; ++Xi[1]) 
                this->Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) 
                    = Flux(Xi[1]-1) - Flux(Xi[1]) + Dj(Xi[1]);

          }
    }
//    CheckDensity(Distribution,"YEnd");
    
    this->boundary->exchangeY(this->Distribution);
    this->boundary->exchangeX(this->Distribution);
//    CheckDensity(Distribution,"YDone");
    
}

template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_x(double timestep) {
        
//    CheckDensity(Distribution,"VXStart");
    errmsg = false;
    
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[2];
    int bvx = UBound[2];
    
    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    VelocityD Vflux;
    
    int j_old; 
    for (Xi[0] = LBound[0]+2; Xi[0] < UBound[0]-1; ++Xi[0]) 
      for (Xi[1] = LBound[1]+2; Xi[1] < UBound[1]-1; ++Xi[1]) {
             
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            Vflux[1] = Vi[1];
            Vflux[2] = Vi[2];
            j_old = lvx-1;
             
            for (Vi[0] = lvx; Vi[0] < bvx; ++Vi[0]) {
                Vflux[0] = Vi[0]+0.5;
                VelocityD Vel = this->velocity(Vflux);
                double F = this->ForceX(Xi,Vel,-timestep);
                double deltavx = F/this->deltaVx();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                int j = max(min(Vi[0]+deltaI, bvx+1), lvx-1);
                // go deltaI right 

//                if ( ((deltaI<0) || (deltaI>1)) && !errmsg ) 
//                {
//                  cerr << "Vx: deltaI out of bounds:" << deltaI << endl;
//                  errmsg = true;
//                }

                double &djref = Dj(Vi[0]);
                
                djref = 0;
                for (int jj=j_old+1; jj<=min(j,bvx) ; ++jj)
                  djref += this->Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);
                  
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
              Dj(bvx) += this->Distribution(Xi[0], Xi[1], jj, Vi[1], Vi[2]);

            
            double &dl = this->Distribution(Xi[0], Xi[1], lvx, Vi[1], Vi[2]);
            Flux(lvx) = min(Flux(lvx), Dj(lvx));
            dl  = - Flux(lvx) + Dj(lvx);
    
            
            for (Vi[0] = lvx+1; Vi[0] < bvx; ++Vi[0]) {
                double &dd = this->Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                dd = Flux(Vi[0]-1) - Flux(Vi[0]) + Dj(Vi[0]);
            }
            
            double &dh = this->Distribution(Xi[0], Xi[1], bvx, Vi[1], Vi[2]);
            Flux(bvx-1) = max(Flux(bvx-1), -Dj(bvx));
            dh  = Flux(bvx-1) + Dj(bvx);
              
        }
    }
//    CheckDensity(Distribution,"VXEnd");
    this->boundary->exchangeX(this->Distribution);
    this->boundary->exchangeY(this->Distribution);
//    CheckDensity(Distribution,"VXDone");
}

template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_y(double timestep) {
        
//    CheckDensity(Distribution,"VYStart");
    errmsg = false;
    
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[3];
    int bvx = UBound[3];

    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    VelocityD Vflux;

    int j_old; 
    
    for (Xi[0] = LBound[0]+2; Xi[0] < UBound[0]-1; ++Xi[0]) 
      for (Xi[1] = LBound[1]+2; Xi[1] < UBound[1]-1; ++Xi[1]) {
               
        for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
          for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            Vflux[0] = Vi[0];
            Vflux[2] = Vi[2];
            
            j_old = lvx-1;
            
            for (Vi[1] = lvx; Vi[1] < bvx; ++Vi[1]) {
                Vflux[1] = Vi[1]+0.5;
                VelocityD Vel = this->velocity(Vflux);
                
                double F = this->ForceY(Xi,Vel,-timestep);
                double deltavx = F/this->deltaVy();
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

                double &djref = Dj(Vi[1]);
                
                djref = 0;
                
                for (int jj=j_old+1; jj<=min(j,bvx); ++jj)
                  djref += this->Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);
                  
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
              Dj(bvx) += this->Distribution(Xi[0], Xi[1], Vi[0], jj, Vi[2]);
  
            double &dl = this->Distribution(Xi[0], Xi[1], Vi[0], lvx, Vi[2]);
            Flux(lvx) = min(Flux(lvx), Dj(lvx));
            dl  = - Flux(lvx) + Dj(lvx);
            
            for (Vi[1] = lvx+1; Vi[1] < bvx; ++Vi[1]) {
                double &dd = this->Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                dd = Flux(Vi[1]-1) - Flux(Vi[1]) + Dj(Vi[1]);
            }

            double &dh = this->Distribution(Xi[0], Xi[1], Vi[0], bvx, Vi[2]);
            Flux(bvx-1) = max(Flux(bvx-1), -Dj(bvx));
            dh = Flux(bvx-1) + Dj(bvx);
        }
    }
//    CheckDensity(Distribution,"VYEnd");
    
    this->boundary->exchangeX(this->Distribution);
    this->boundary->exchangeY(this->Distribution);
//    CheckDensity(Distribution,"VYDone");
}


template<class ForceField>
void PosFluxCons3rdOrder<ForceField>
        ::advanceVel_z(double timestep) {
        
//    CheckDensity(Distribution,"VZStart");
    errmsg = false;
    
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();

    PositionI Xi;
    VelocityI Vi;
    
    int lvx = LBound[4];
    int bvx = UBound[4];

    NumMatrix<double, 1> Flux(&lvx,&bvx);
    NumMatrix<double, 1> Dj(&lvx,&bvx);
    
    VelocityD Vflux;

    int j_old; // check initialization inside the loop!!
    
//    cerr << "Advancing in vz\n";
    for (Xi[0] = LBound[0]+2; Xi[0] < UBound[0]-1; ++Xi[0]) 
      for (Xi[1] = LBound[1]+2; Xi[1] < UBound[1]-1; ++Xi[1]) {
      
        for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
          for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) {
            Vflux[0] = Vi[0];
            Vflux[1] = Vi[1];
            
            j_old = lvx-1;
             
            for (Vi[2] = lvx; Vi[2] < bvx; ++Vi[2]) {
            
                Vflux[2] = Vi[2]+0.5;
                VelocityD Vel = this->velocity(Vflux);
                
                double F = this->ForceZ(Xi,Vel,-timestep);
                double deltavx = F/this->deltaVz();
                int deltaI = int(floor(deltavx)+1);  
                double alpha = deltavx - deltaI + 1;

                // go deltaI right 
                int j = max(min(Vi[2]+deltaI, bvx+1), lvx-1);
                
//                if ( ((deltaI<0) || (deltaI>1)) && !errmsg ) 
//                {
//                  cerr << "Vz: deltaI out of bounds:" << deltaI << endl;
//                  errmsg = true;
//                }
                
                double &djref = Dj(Vi[2]);
                
                djref = 0;
                int jjbound = min(j,bvx);
                for (int jj=j_old+1; jj<=jjbound; ++jj)
                  djref += this->Distribution(Xi[0], Xi[1],Vi[0], Vi[1], jj);

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
              Dj(bvx) += this->Distribution(Xi[0], Xi[1], Vi[0], Vi[1], jj);
  

            double &dl = this->Distribution(Xi[0], Xi[1], Vi[0], Vi[1], lvx);
            
            Flux(lvx) = min(Flux(lvx), Dj(lvx));
            
            dl = - Flux(lvx) + Dj(lvx);
            
            for (Vi[2] = lvx+1; Vi[2] < bvx; ++Vi[2]) {
                double &dd = this->Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
                dd = Flux(Vi[2]-1) - Flux(Vi[2]) + Dj(Vi[2]);
            }
            double &dh = this->Distribution(Xi[0], Xi[1], Vi[0], Vi[1], bvx);
            
            Flux(bvx-1) = max(Flux(bvx-1), -Dj(bvx));
            
            dh = Flux(bvx-1) + Dj(bvx);
        }
    }
//    CheckDensity(Distribution,"VZEnd");
    
    this->boundary->exchangeX(this->Distribution);
    this->boundary->exchangeY(this->Distribution);
//    CheckDensity(Distribution,"VZDone");

}

//#define FULL_LIMIT
//#define NO_LIMIT

template<class ForceField>
inline double PosFluxCons3rdOrder<ForceField>
  ::epsilonLeft
    (
      double fj, 
      double fjp
    ) 
{

#ifdef NO_LIMIT
  return 1.;
#endif

  double fdiff = fjp-fj;
  double fdiff_safe = fdiff + double(fdiff==0);

  double result[3];
  result[0] = 1.;
  result[1] = 2*fj/fdiff_safe;

  int res1 = (2*fj<fdiff);

#ifdef FULL_LIMIT
  double fexc=2*(f_infty-fj);    
  result[2] = -fexc/fdiff_safe;
  int res2 = (fexc < (-fdiff) );
  return result[res1*(1-res2) + 2*(1-res1)*res2];
#else     
  return result[res1];
#endif
}

template<class ForceField>
inline double PosFluxCons3rdOrder<ForceField>
  ::epsilonRight
    (
      double fj, 
      double fjm
    ) 
{

#ifdef NO_LIMIT
  return 1.;
#endif

  double fdiff = fjm-fj;
  double fdiff_safe = fdiff + double(fdiff==0);

  double result[3];
  result[0] = 1.;
  result[1] = 2*fj/fdiff_safe;
  int res1 = (2*fj<fdiff);

#ifdef FULL_LIMIT
  double fexc=2*(f_infty-fj); 
  result[2] = -fexc/fdiff_safe;
  int res2 = (fexc < (-fdiff) );
  return result[res1*(1-res2) + 2*(1-res1)*res2];
#else
  return result[res1];
#endif
}


template<class ForceField>
double PosFluxCons3rdOrder<ForceField>::interpolateX(
                const PositionI &Xi,
                const VelocityI &Vi,
                int j,
                double alpha) {
    double fjp = this->Distribution(j+1,Xi[1],Vi[0],Vi[1],Vi[2]);
    double fj  = this->Distribution(j  ,Xi[1],Vi[0],Vi[1],Vi[2]);
    double fjm = this->Distribution(j-1,Xi[1],Vi[0],Vi[1],Vi[2]);
    
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
    double fjp = this->Distribution(Xi[0], j+1,Vi[0],Vi[1],Vi[2]);
    double fj  = this->Distribution(Xi[0], j  ,Vi[0],Vi[1],Vi[2]);
    double fjm = this->Distribution(Xi[0], j-1,Vi[0],Vi[1],Vi[2]);
    
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
    double fj  = this->Distribution(Xi[0],Xi[1],j  ,Vi[1],Vi[2]);
    double fjp, fjm;
    
//     double fjparr[3], fjmarr[3];
//     
//     fjp = Distribution(Xi[0],Xi[1],j+1,Vi[1],Vi[2]);
//     fjm = Distribution(Xi[0],Xi[1],j-1,Vi[1],Vi[2]);
//     
//     fjparr[0] = fjp;
//     fjparr[1] = 0.0;
//     fjparr[2] = 2*fj-fjm;
//     fjmarr[0] = fjm;
//     fjmarr[1] = 0.0;
//     fjmarr[2] = 2*fj-fjp;
//     
//     int logicp1 = int(j >= Distribution.getHigh()[2]);
//     int logicp2 = int( fjparr[2] > 0 );
//     int logicm1 = int(j <= Distribution.getLow()[2]);
//     int logicm2 = int( fjmarr[2] > 0 );
//     
//     fjp = fjparr[ logicp1*(1+logicp2) ];
//     fjm = fjmarr[ logicm1*(1+logicm2) ];
  
  if (j<=this->Distribution.getLow()[2])
  {
    fjp = this->Distribution(Xi[0],Xi[1],j+1,Vi[1],Vi[2]);
    fjm = max(0.0,2*fj-fjp);
  } else
  if (j>=this->Distribution.getHigh()[2])
  {
    fjm = this->Distribution(Xi[0],Xi[1],j-1,Vi[1],Vi[2]);
    fjp = max(0.0,2*fj-fjm);
  } else
  {
    fjp = this->Distribution(Xi[0],Xi[1],j+1,Vi[1],Vi[2]);
    fjm = this->Distribution(Xi[0],Xi[1],j-1,Vi[1],Vi[2]);
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
    double fj  = this->Distribution(Xi[0],Xi[1],Vi[0],j  ,Vi[2]);
    double fjp, fjm;
    if (j<=this->Distribution.getLow()[3])
    {
      fjp = this->Distribution(Xi[0],Xi[1],Vi[0],j+1,Vi[2]);
      fjm = max(0.0,2*fj-fjp);
    } else
    if (j>=this->Distribution.getHigh()[3])
    {
      fjm = this->Distribution(Xi[0],Xi[1],Vi[0],j-1,Vi[2]);
      fjp = max(0.0,2*fj-fjm);
    } else
    {
      fjp = this->Distribution(Xi[0],Xi[1],Vi[0],j+1,Vi[2]);
      fjm = this->Distribution(Xi[0],Xi[1],Vi[0],j-1,Vi[2]);
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
    double fj  = this->Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j  );
    double fjp, fjm;
    if (j<=this->Distribution.getLow()[4])
    {
      fjp = this->Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j+1);
      fjm = max(0.0,2*fj-fjp);
    } else
    if (j>=this->Distribution.getHigh()[4])
    {
      fjm = this->Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j-1);
      fjp = max(0.0,2*fj-fjm);
    } else
    {
      fjp = this->Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j+1);
      fjm = this->Distribution(Xi[0],Xi[1],Vi[0],Vi[1],j-1);
    }

    double epsl = epsilonLeft(fj,fjp); // Slope Corrector
    double epsr = epsilonRight(fj,fjm); // Slope Corrector
    return (1-alpha)*(fj
                +(1/6.)*alpha*(alpha + 1)*epsl*(-fj+fjp)
                -(1/6.)*alpha*(alpha - 2)*epsr*(-fjm+fj));
}

