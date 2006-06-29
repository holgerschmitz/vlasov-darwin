
BoundMatrix::BoundMatrix(VlasovDist &dist_, DirectionType direction_)
  dist(dist_),
  direction(direction_)
{
  setDirection(direction);
}

void BoundMatrix::setDirection(DirectionType direction_)
{
  direction = direction_;
  int *dlo = dist.getLow();
  int *dhi = dist.getHigh();
  
  for (int d = 0; d < 5; d++) {
    lo[d] = dlo[d];
    hi[d] = dhi[d]; 
    dims[d] = dhi[d] - dlo[d] + 1;
  }
  
  switch (direction)
  {
    case left:
      hi[0] = lo[0]+1;
      dims[0] = 2;
      break;
    case right:
      lo[0] = hi[0]-1;
      dims[0] = 2;
      break;
    case bottom:
      hi[1] = lo[1]+1;
      dims[1] = 2;
      break;
    case top:
      lo[1] = hi[1]-1;
      dims[1] = 2;
      break;
  }
}

inline
void BoundMatrix::copyVelocity(int fromi, int fromj, int toi, int toj)
{
  for (int vx=lo[2]; vx<=hi[2]; ++vx)
    for (int vy=lo[3]; vy<=hi[3]; ++vy)
      for (int vz=lo[4]; vz<=hi[4]; ++vz)
      {
        dist(toi,toj,vx,vy,vz) = dist(fromi,fromj,vx,vy,vz);
      }
}

void BoundMatrix::copyBoundary()
{
  switch (direction)
  {
    case left:
      int i0 = lo[0];
      int i1 = lo[0]+1;
      int i2 = lo[0]+2;
      int i3 = lo[0]+3;
      for (int j=lo[1]; j<=hi[1]; ++j)
      {
        copyVelocity(i3,j,i0,j);
        copyVelocity(i2,j,i1,j);
      }
      break;
      
    case right:
      int i0 = hi[0];
      int i1 = hi[0]-1;
      int i2 = hi[0]-2;
      int i3 = hi[0]-3;
      for (int j=lo[1]; j<=hi[1]; ++j)
      {
        copyVelocity(i3,j,i0,j);
        copyVelocity(i2,j,i1,j);
      }
      break;
      
    case bottom:
      int j0 = lo[1];
      int j1 = lo[1]+1;
      int j2 = lo[1]+2;
      int j3 = lo[1]+3;
      for (int i=lo[0]; i<=hi[0]; ++i)
      {
        copyVelocity(i,j3,i,j0);
        copyVelocity(i,j2,i,j1);
      }
      break;

    case top:
      int j0 = hi[1];
      int j1 = hi[1]-1;
      int j2 = hi[1]-2;
      int j3 = hi[1]-3;
      for (int i=lo[0]; i<=hi[0]; ++i)
      {
        copyVelocity(i,j3,i,j0);
        copyVelocity(i,j2,i,j1);
      }
      break;
  }
}


OpenBoundForce::OpenBoundForce
  (
    ForceFieldBase &base_
  )
  : base(base_),
    Distribution(base.Distribution)
{
  dx = base.dx;
  dv[0] = base.deltaVx();
  dv[1] = base.deltaVy();
  dv[2] = base.deltaVz();
  
  // only the first two indices are used for resizing
  veloffset.resize(Distribution.getLow(), Distribution.getHigh);
}

void OpenBoundForce::setDirection(BoundMatrix::DirectionType direction)
{
  Distribution.setDirection(direction);

  // only the first two indices are used for resizing
  veloffset.resize(Distribution.getLow(), Distribution.getHigh);
}

VelocityD OpenBoundForce::Force(const PositionI &Pos, 
                const VelocityD &Vel,
                double dt)
{
  return veloffset(i,j);
}

double OpenBoundForce::ForceX(const PositionI &Pos, 
              const VelocityD &Vel,
              double dt)
{
  return veloffset(i,j)[0];
}

double OpenBoundForce::ForceY(const PositionI &Pos, 
              const VelocityD &Vel,
              double dt)
{
  return veloffset(i,j)[1];
}

double OpenBoundForce::ForceZ(const PositionI &Pos, 
              const VelocityD &Vel,
              double dt)
{
  return veloffset(i,j)[2];
}


void OpenBoundForce::setMoments(double density_, VelocityD current_)
{
  density = density_;
  current = current_;
}

void OpenBoundForce::initialize()
{
  int *lo = Distribution.getLow();
  int *hi = Distribution.getHigh();

  copyBoundary();
  for (int i=lo[0]; i<=hi[0]; ++i)
    for (int j=lo[1]; j<=hi[1]; ++j)
    {
      double realden = integrateDensity(i,j);
      multiplyDensity(i, j, density/realden);
      VelocityD realcurr = integrateCurrent(i,j);
      veloffset(i,j) = (current-realcurr)/density;
      veloffset(i,j) *= 2;
    }
}

double OpenBoundForce::integrateDensity(int i, int j)
{
  int *lo = Distribution.getLow();
  int *hi = Distribution.getHigh();
  
  double den = 0;
  for (int k=lo[2]; k<=hi[2]; ++k)
    for (int l=lo[3]; l<=hi[3]; ++l)
      for (int m=lo[4]; m<=hi[4]; ++m) 
        den += Distribution(i,j,k,l,m);
  return den;
}

void OpenBoundForce::multiplyDist(int i, int j, double factor)
{
  int *lo = Distribution.getLow();
  int *hi = Distribution.getHigh();
  
  double den = 0;
  for (int k=lo[2]; k<=hi[2]; ++k)
    for (int l=lo[3]; l<=hi[3]; ++l)
      for (int m=lo[4]; m<=hi[4]; ++m) 
        Distribution(i,j,k,l,m) *= factor;
  return den;
}

VelocityD OpenBoundForce::integrateCurrent(int i, int j) 
{
  int *lo = Distribution.getLow();
  int *hi = Distribution.getHigh();

  double jxval=0, jyval=0, jzval=0;
  VelocityI vi;
  VelocityD V;
  
  for (int vi[0]=lo[2]; vi[0]<=hi[2]; ++vi[0])
    for (int vi[1]=lo[3]; vi[1]<=hi[3]; ++vi[1])
      for (int vi[2]=lo[4]; vi[2]<=hi[4]; ++vi[2])
      {
        V = base.velocity(vi);
        d = Distribution(i,j,vi[0],vi[1],vi[2]);
        
        jxval += V[0]*d;
        jyval += V[1]*d;
        jzval += V[2]*d;
      }
  
  return VelocityD(jxval,jyval, jzval);
}


OpenBound::OpenBound(ForceFieldBase& base, BoundMatrix::DirectionType direction)
  : PosFluxCons3rdOrder(base)
{
  setDirection(direction);
}

void OpenBound::apply()
{
  this->initialize();
  
  // the time step does not play a role
  advanceVel_x(1.0);
  advanceVel_y(1.0);
  advanceVel_z(1.0);
}
