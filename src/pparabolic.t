
#include "numeric.h"

template<class ForceField>
void PParabolicScheme<ForceField>
        ::advanceSpace_x(double timestep) 
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  PositionI Xi;
  VelocityI Vi;
  
  int lx = LBound[0];
  int bx = UBound[0];

  du.resize(&lx,&bx);
  uh.resize(&lx,&bx);
  uL.resize(&lx,&bx);
  uR.resize(&lx,&bx);
  u6.resize(&lx,&bx);
  Phi.resize(&lx,&bx);
  
  NumMatrix<double, 1> v(&lx,&bx);
  NumMatrix<double, 1> f(&lx,&bx);
  BoundX bound(boundary);
  
  for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) {
    double deltax = -timestep*velocity(Vi)[0]/dx[0];

    for (Xi[1] = LBound[1]; Xi[1] <= UBound[1]; ++Xi[1])
      for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
        for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            
             
          // Copy into temporary field
          for (Xi[0] = lx; Xi[0] <= bx; ++Xi[0]) 
          {
            f(Xi[0]) = Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]);
            v(Xi[0]) = deltax;
          }
          
          oneDimPpm(f,v,lx,bx, &bound);
          
          // Copy back from temporary field
          for (Xi[0] = lx; Xi[0] <= bx; ++Xi[0]) 
          {
            Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) = f(Xi[0]);
          }
        }
  }
            
}

template<class ForceField>
void PParabolicScheme<ForceField>
        ::advanceSpace_y(double timestep) 
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  PositionI Xi;
  VelocityI Vi;
  
  int lx = LBound[1];
  int bx = UBound[1];

  du.resize(&lx,&bx);
  uh.resize(&lx,&bx);
  uL.resize(&lx,&bx);
  uR.resize(&lx,&bx);
  u6.resize(&lx,&bx);
  Phi.resize(&lx,&bx);
  
  NumMatrix<double, 1> v(&lx,&bx);
  NumMatrix<double, 1> f(&lx,&bx);
  BoundY bound(boundary);
  
  for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) {
    double deltax = -timestep*velocity(Vi)[1]/dx[1];

    for (Xi[0] = LBound[0]; Xi[0] <= UBound[0]; ++Xi[0])
      for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
        for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            
             
          // Copy into temporary field
          for (Xi[1] = lx; Xi[1] <= bx; ++Xi[1]) 
          {
            f(Xi[1]) = Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]);
            v(Xi[1]) = deltax;
          }
          
          oneDimPpm(f,v,lx,bx, &bound);
          
          // Copy back from temporary field
          for (Xi[1] = lx; Xi[1] <= bx; ++Xi[1]) 
          {
            Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) = f(Xi[1]);
          }
        }
  }
            
}

template<class ForceField>
void PParabolicScheme<ForceField>
        ::advanceVel_x(double timestep) 
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  PositionI Xi;
  VelocityI Vi;
  
  int lvx = LBound[2];
  int bvx = UBound[2];

  du.resize(&lvx,&bvx);
  uh.resize(&lvx,&bvx);
  uL.resize(&lvx,&bvx);
  uR.resize(&lvx,&bvx);
  u6.resize(&lvx,&bvx);
  Phi.resize(&lvx,&bvx);
  
  NumMatrix<double, 1> v(&lvx,&bvx);
  NumMatrix<double, 1> f(&lvx,&bvx);

  VelocityD Vflux;

  BoundV bound;
  
  for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
    for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1])
           
      for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) 
        for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {

          Vflux[1] = Vi[1];
          Vflux[2] = Vi[2];
            
          for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0]) {
            
            Vflux[0] = Vi[0]+0.5;
            VelocityD Vel = velocity(Vflux);
            double F = ForceX(Xi,Vel,timestep);
             
            f(Vi[0]) = Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]);
            v(Vi[0]) = -F/deltaVx();
          }            
          
          oneDimPpm(f,v,lvx,bvx, &bound);
            
          // Copy back from temporary field
          for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0]) 
          {
            Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) = f(Vi[0]);
          }
        }
           
}

template<class ForceField>
void PParabolicScheme<ForceField>
        ::advanceVel_y(double timestep) 
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  PositionI Xi;
  VelocityI Vi;
  
  int lvx = LBound[3];
  int bvx = UBound[3];

  du.resize(&lvx,&bvx);
  uh.resize(&lvx,&bvx);
  uL.resize(&lvx,&bvx);
  uR.resize(&lvx,&bvx);
  u6.resize(&lvx,&bvx);
  Phi.resize(&lvx,&bvx);
  
  NumMatrix<double, 1> v(&lvx,&bvx);
  NumMatrix<double, 1> f(&lvx,&bvx);

  VelocityD Vflux;

  BoundV bound;
  
  for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
    for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1])
           
      for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
        for (Vi[2] = LBound[4]; Vi[2] <= UBound[4]; ++Vi[2]) {
            
          Vflux[0] = Vi[0];
          Vflux[2] = Vi[2];

          for (Vi[1] = lvx; Vi[1] <= bvx; ++Vi[1]) {
            
            Vflux[1] = Vi[1]+0.5;
            VelocityD Vel = velocity(Vflux);
            double F = ForceY(Xi,Vel,timestep);
             
            f(Vi[1]) = Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]);
            v(Vi[1]) = -F/deltaVx();
          }            
          
          oneDimPpm(f,v,lvx,bvx, &bound);
            
          // Copy back from temporary field
          for (Vi[1] = lvx; Vi[1] <= bvx; ++Vi[1]) 
          {
            Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) = f(Vi[1]);
          }
        }
           
}

template<class ForceField>
void PParabolicScheme<ForceField>
        ::advanceVel_z(double timestep) 
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  PositionI Xi;
  VelocityI Vi;
  
  int lvx = LBound[4];
  int bvx = UBound[4];

  du.resize(&lvx,&bvx);
  uh.resize(&lvx,&bvx);
  uL.resize(&lvx,&bvx);
  uR.resize(&lvx,&bvx);
  u6.resize(&lvx,&bvx);
  Phi.resize(&lvx,&bvx);
  
  NumMatrix<double, 1> v(&lvx,&bvx);
  NumMatrix<double, 1> f(&lvx,&bvx);

  VelocityD Vflux;

  BoundV bound;
  
  for (Xi[0] = LBound[0]+1; Xi[0] < UBound[0]; ++Xi[0]) 
    for (Xi[1] = LBound[1]+1; Xi[1] < UBound[1]; ++Xi[1])
           
      for (Vi[0] = LBound[2]; Vi[0] <= UBound[2]; ++Vi[0]) 
        for (Vi[1] = LBound[3]; Vi[1] <= UBound[3]; ++Vi[1]) {
            
          Vflux[0] = Vi[0];
          Vflux[1] = Vi[1];

          for (Vi[2] = lvx; Vi[2] <= bvx; ++Vi[2]) {
            
            Vflux[2] = Vi[2]+0.5;
            VelocityD Vel = velocity(Vflux);
            double F = ForceZ(Xi,Vel,timestep);
             
            f(Vi[2]) = Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]);
            v(Vi[2]) = -F/deltaVx();
          }            
          
          oneDimPpm(f,v,lvx,bvx, &bound);
            
          // Copy back from temporary field
          for (Vi[2] = lvx; Vi[2] <= bvx; ++Vi[2]) 
          {
            Distribution(Xi[0] ,Xi[1] , Vi[0], Vi[1], Vi[2]) = f(Vi[2]);
          }
        }
           
}

template<class ForceField>
void PParabolicScheme<ForceField>
        ::oneDimPpm(NumMatrix<double, 1> &u,
                    NumMatrix<double, 1> &v,
                    int lx, int hx,
                    BoundaryBase *bound)
{


  int lxp = lx+1, lxP = lx+2;
  int hxm = hx-1;

  for(int m = lxp; m < hx; ++m)
  {
    if( (u(m+1)-u(m))*(u(m)-u(m-1)) > 0.0)
    {
      double diff = 0.5 * (u(m+1) - u(m-1));
      du(m) = sgn(diff)*min(
                fabs(diff),
                2.0*min( fabs(u(m)-u(m-1)), fabs(u(m+1)-u(m)) )
              );
    }
    else
      du(m) = 0.0;
  }
  
  bound->apply(du);

  for(int m = lx; m < hx; m++)
    uh(m) = (u(m+1)+u(m))/2.0-(du(m+1)-du(m))/6.0;

  for(int m = lxp; m < hx; m++)
  {
    uL(m) = uh(m - 1);
    uR(m) = uh(m);
  }

  for(int m = lxp; m < hx; m++)
  {
    du(m) = uR(m) - uL(m);
    u6(m) = 6.0*(u(m) - 0.5 * (uL(m) + uR(m)));
    
    if((uR(m)-u(m))*(u(m)-uL(m)) <= 0.0)
    {
      uR(m) = u(m);
      uL(m) = u(m);
    }
    else if( du(m) * u6(m) > du(m) * du(m)) 
      uL(m) = 3.0 * u(m) - 2.0 * uR(m);
    else if(-du(m) * du(m) > du(m) * u6(m)) 
      uR(m) = 3.0 * u(m) - 2.0 * uL(m);
      
    du(m) = uR(m) - uL(m);
    u6(m) = 6.0*(u(m) - 0.5 * (uL(m) + uR(m)));
  }

  for(int m = lxp; m < hxm; m++)
  {
    if(v(m) >= 0.0) 
    {
      Phi(m) = v(m)*(uR(m)-0.5*v(m)*(du(m)-(1.0-2.0/3.0*v(m))*u6(m)));
    }
    else 
    {
      Phi(m) = v(m)*(uL(m+1)-0.5*v(m)*(du(m+1)+(1.0+2.0/3.0*v(m))*u6(m+1)));
    }
  }
  
  for(int m = lxP; m < hxm; m++) 
    u(m) -= (Phi(m) - Phi(m-1));

  bound->apply(u);
}

//int ppm(double v_mu, long int N, int bc_flag, double* ppm_buffer)
//{
/*  Piecewise Parabolic Method for dz/dt + v dz/dx = 0
   z_in[0..N] = array of pointers to input data
   z_out[0..N] = array of pointers to output data
   mu = dt/dx
   v = advection velocity
   v_mu = v_mu
    N = number of data points
    bc_flag = 0 if and only if data has periodic boundary conditions */

   /* double ppm_buffer[7*N+m_ext+m_ext]; */

// long int m; /* local dummy variable */
// const long int S = N + m_ext + m_ext;
//  double *u, *uh, *uL, *uR, *du, *u6, *Phi;
//  m = 0;
// u = &(ppm_buffer[m]);
//  m += S;
//  uh = &(ppm_buffer[m]);
// m += S;
//  uL = &(ppm_buffer[m]);
// m += S;
//  uR = &(ppm_buffer[m]);
// m += S;
//  du = &(ppm_buffer[m]);
// m += S;
//  u6 = &(ppm_buffer[m]);
// m += S;
//  Phi = &(ppm_buffer[m]);

// if(fabs(v_mu) > 1.0) pre_shift(&v_mu, S, bc_flag, ppm_buffer);
//
//  for(m = -m_ext+1; m < N+m_ext-1; m++)
//  {
//    if((u[m+m_ext+1]-u[m+m_ext])*(u[m+m_ext]-u[m+m_ext-1]) > 0.0)
//    {
//      du[m+m_ext] = 0.5 * (u[m+m_ext+1] - u[m+m_ext-1]);
//     du[m+m_ext] = sgn(du[m+m_ext])*min(fabs(du[m+m_ext]),2.0*min(fabs(u[m+m_ext]-u[m+m_ext-1]),fabs(u[m+m_ext+1]-u[m+m_ext])));
//    }
//    else
//     du[m+m_ext] = 0.0;
// }

//  if(bc_flag == 0)
// {
//    du[0] = du[N];
//   du[S - 1] = du[m_ext + m_ext - 1];
// }
//  else
// {
//    du[0] = 0.0;
//   du[S - 1] = 0.0;
// }

// for(m = -m_ext; m < N+m_ext-1; m++)
//    uh[m+m_ext] = (u[m+m_ext+1]+u[m+m_ext])/2.0-(du[m+m_ext+1]-du[m+m_ext])/6.0;
//
//  if(bc_flag == 0)
//   uh[S-1] = uh[m_ext + m_ext-1];
// else
//   uh[S - 1] = 0.0;

//  for(m = -m_ext+1; m < N+m_ext; m++)
//  {
//    uL[m + m_ext] = uh[m + m_ext - 1];
//   uR[m + m_ext] = uh[m + m_ext];
// }
//  if(bc_flag == 0)
// {
//    uL[0] = uL[N];
//   uR[0] = uR[N];
// }
//  else
// {
//    uL[0] = 0.0;
//   uR[0] = 0.0;
// }

// for(m = -m_ext; m < N+m_ext; m++)
//  {
//    du[m+m_ext] = uR[m+m_ext] - uL[m+m_ext];
//   u6[m+m_ext] = 6.0*(u[m+m_ext] - 0.5 * (uL[m+m_ext] + uR[m+m_ext]));
//    if((uR[m+m_ext]-u[m+m_ext])*(u[m+m_ext]-uL[m+m_ext]) <= 0.0)
//   {
//      uR[m+m_ext] = u[m+m_ext];
//      uL[m+m_ext] = u[m+m_ext];
//    }
//    else if( du[m+m_ext] * u6[m+m_ext] > du[m+m_ext] * du[m+m_ext]) uL[m+m_ext] = 3.0 * u[m+m_ext] - 2.0 * uR[m+m_ext];
//    else if(-du[m+m_ext] * du[m+m_ext] > du[m+m_ext] * u6[m+m_ext]) uR[m+m_ext] = 3.0 * u[m+m_ext] - 2.0 * uL[m+m_ext];
//    du[m+m_ext] = uR[m+m_ext] - uL[m+m_ext];
//   u6[m+m_ext] = 6.0*(u[m+m_ext] - 0.5 * (uL[m+m_ext] + uR[m+m_ext]));
//  }

// if(v_mu >= 0.0) for(m = -1; m < N; m++)
//  {
//    Phi[m+m_ext] = v_mu*(uR[m+m_ext]-0.5*v_mu*(du[m+m_ext]-(1.0-2.0/3.0*v_mu)*u6[m+m_ext]));
// }
//  else for(m = -1; m < N; m++)
// {
//    Phi[m+m_ext] = v_mu*(uL[m+m_ext+1]-0.5*v_mu*(du[m+m_ext+1]+(1.0+2.0/3.0*v_mu)*u6[m+m_ext+1]));
// }

// for(m = 0; m < N; m++) u[m+m_ext] -= (Phi[m+m_ext] - Phi[m+m_ext-1]);
//
// return 0;
//}
