//================================================================
//========       FiniteDiffScheme        =========================
//================================================================

template<class ForceField>
void FiniteDiffScheme<ForceField>::initDist()
{
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
    
  TempDist.resize(LBound,UBound);
  CurrentDist.resize(LBound,UBound);
}

template<class ForceField>
void FiniteDiffScheme<ForceField>
        ::advanceAll(double timestep) {
    
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();

                                              //  T | C | D |
    this->Distribution.swap(CurrentDist);     //  T | D | C |
    this->Distribution.swap(TempDist);        //  D | T | C |
                                              // after the following we have   
                                              //    | T | C | D
  
    PositionI Xi;
    VelocityI Vi;
    
    int lx = LBound[0]+1;
    int bx = UBound[0]-1;
    int ly = LBound[1]+1;
    int by = UBound[1]-1;
    int lvx = LBound[2]+1;
    int bvx = UBound[2]-1;
    int lvy = LBound[3]+1;
    int bvy = UBound[3]-1;
    int lvz = LBound[4]+1;
    int bvz = UBound[4]-1;
    
    double dvx = this->deltaVx();
    double dvy = this->deltaVy();
    double dvz = this->deltaVz();
    
    for (Vi[2] = lvz; Vi[2] <= bvz; ++Vi[2])
      for (Vi[1] = lvy; Vi[1] <= bvy; ++Vi[1])	
        for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0])
          for (Xi[1] = ly;  Xi[1] <= by;  ++Xi[1])
            for (Xi[0] = lx;  Xi[0] <= bx;  ++Xi[0])
            {
              VelocityD Vel = this->velocity(Vi);
              VelocityD F = this->Force(Xi,Vel,timestep);
              
              double dfx = (CurrentDist(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          - CurrentDist(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])) / (2*this->dx[0]);

              double dfy = (CurrentDist(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          - CurrentDist(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])) / (2*this->dx[1]);

              double dfvx = (CurrentDist(Xi[0], Xi[1], Vi[0]+1, Vi[1], Vi[2])
                           - CurrentDist(Xi[0], Xi[1], Vi[0]-1, Vi[1], Vi[2])) / (2*dvx);

              double dfvy = (CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1]+1, Vi[2])
                           - CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1]-1, Vi[2])) / (2*dvy);

              double dfvz = (CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+1)
                           - CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-1)) / (2*dvz);
                           
              this->Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) 
                = TempDist(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                   - 2*( timestep*(Vel[0]*dfx + Vel[1]*dfy)
                          + (F[0]*dfvx + F[1]*dfvy + F[2]*dfvz)
                        );
            }
    this->boundary->exchangeX(this->Distribution);
    this->boundary->exchangeY(this->Distribution);
    
}

//===================================================================
//==============FiniteDiffSchemeRK2 =================================
//===================================================================

template<class ForceField>
void FiniteDiffSchemeRK2<ForceField>::initDist()
{
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
    
  D1.resize(LBound,UBound);
}

template<class ForceField>
void FiniteDiffSchemeRK2<ForceField>
        ::forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep)
{
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();
    
    PositionI Xi;
    VelocityI Vi;
    
    int lx = LBound[0]+1;
    int bx = UBound[0]-1;
    int ly = LBound[1]+1;
    int by = UBound[1]-1;
    int lvx = LBound[2]+1;
    int bvx = UBound[2]-1;
    int lvy = LBound[3]+1;
    int bvy = UBound[3]-1;
    int lvz = LBound[4]+1;
    int bvz = UBound[4]-1;
    
    double dvx = this->deltaVx();
    double dvy = this->deltaVy();
    double dvz = this->deltaVz();
    
    for (Vi[2] = lvz; Vi[2] <= bvz; ++Vi[2])
      for (Vi[1] = lvy; Vi[1] <= bvy; ++Vi[1])	
        for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0])
          for (Xi[1] = ly;  Xi[1] <= by;  ++Xi[1])
            for (Xi[0] = lx;  Xi[0] <= bx;  ++Xi[0])
            {
              VelocityD Vel = this->velocity(Vi);
              VelocityD F = this->Force(Xi,Vel,timestep);
              
              double dfx = (rhs(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])) / (2*this->dx[0]);

              double dfy = (rhs(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])) / (2*this->dx[1]);

              double dfvx = (rhs(Xi[0], Xi[1], Vi[0]+1, Vi[1], Vi[2])
                           - rhs(Xi[0], Xi[1], Vi[0]-1, Vi[1], Vi[2])) / (2*dvx);

              double dfvy = (rhs(Xi[0], Xi[1], Vi[0], Vi[1]+1, Vi[2])
                           - rhs(Xi[0], Xi[1], Vi[0], Vi[1]-1, Vi[2])) / (2*dvy);

              double dfvz = (rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+1)
                           - rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-1)) / (2*dvz);
                           
              lhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) 
                += - 2*( timestep*(Vel[0]*dfx + Vel[1]*dfy)
                          + (F[0]*dfvx + F[1]*dfvy + F[2]*dfvz)
                        );
            }
    this->boundary->exchangeX(lhs);
    this->boundary->exchangeY(lhs);
}

template<class ForceField>
void FiniteDiffSchemeRK2<ForceField>
        ::advanceAll(double timestep) 
{   
  switch (RKState)
  {
    case 0:
      forwardDist(D1, this->Distribution, 0.5*timestep);
      RKState = 1;
      break;
    case 1:
      forwardDist(this->Distribution, D1, timestep);
      RKState = 0;
      break;
  }
}




//===================================================================
//==============FiniteDiffSchemeRK4 =================================
//===================================================================

template<class ForceField>
void FiniteDiffSchemeRK3<ForceField>
        ::forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep)
{
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();
    
    PositionI Xi;
    VelocityI Vi;
    
    int lx = LBound[0]+1;
    int bx = UBound[0]-1;
    int ly = LBound[1]+1;
    int by = UBound[1]-1;
    int lvx = LBound[2]+1;
    int bvx = UBound[2]-1;
    int lvy = LBound[3]+1;
    int bvy = UBound[3]-1;
    int lvz = LBound[4]+1;
    int bvz = UBound[4]-1;
    
    double dvx = this->deltaVx();
    double dvy = this->deltaVy();
    double dvz = this->deltaVz();
    
    for (Vi[2] = lvz; Vi[2] <= bvz; ++Vi[2])
      for (Vi[1] = lvy; Vi[1] <= bvy; ++Vi[1])	
        for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0])
          for (Xi[1] = ly;  Xi[1] <= by;  ++Xi[1])
            for (Xi[0] = lx;  Xi[0] <= bx;  ++Xi[0])
            {
              VelocityD Vel = this->velocity(Vi);
              VelocityD F = this->Force(Xi,Vel,timestep);
              
              double dfx = (rhs(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])) / (2*this->dx[0]);

              double dfy = (rhs(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])) / (2*this->dx[1]);

              double dfvx = (rhs(Xi[0], Xi[1], Vi[0]+1, Vi[1], Vi[2])
                           - rhs(Xi[0], Xi[1], Vi[0]-1, Vi[1], Vi[2])) / (2*dvx);

              double dfvy = (rhs(Xi[0], Xi[1], Vi[0], Vi[1]+1, Vi[2])
                           - rhs(Xi[0], Xi[1], Vi[0], Vi[1]-1, Vi[2])) / (2*dvy);

              double dfvz = (rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+1)
                           - rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-1)) / (2*dvz);
                           
              lhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) 
                += - ( timestep*(Vel[0]*dfx + Vel[1]*dfy)
                          + (F[0]*dfvx + F[1]*dfvy + F[2]*dfvz)
                        );
            }
    this->boundary->exchangeX(lhs);
    this->boundary->exchangeY(lhs);
    
    double denErr = this->densityError(lhs)/numDensity;
    this->correctDensityError(denErr, lhs);
}

template<class ForceField>
void FiniteDiffSchemeRK3<ForceField>::initDist()
{
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
    
  T1Dist.resize(LBound,UBound);
  T2Dist.resize(LBound,UBound);

  for (int i=LBound[0]; i<=UBound[0]; ++i)
    for (int j=LBound[1]; j<=UBound[1]; ++j)
      for (int k=LBound[2]; k<=UBound[2]; ++k)
        for (int l=LBound[3]; l<=UBound[3]; ++l) 
          for (int m=LBound[4]; m<=UBound[4]; ++m)
          {
            T1Dist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
          }
  if (this->boundary->master()) std::cout << "Initializing Density Error\n";
  numDensity = this->densityError(this->Distribution);
}


template<class ForceField>
void FiniteDiffSchemeRK3<ForceField>
        ::advanceAll(double timestep) 
{
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();

  switch (RKState) {
      case 0: forwardDist(this->Distribution, T1Dist, timestep/3.);         
              this->Distribution.swap(T1Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                      {
                        T2Dist(i,j,k,l,m) = 0.25*T1Dist(i,j,k,l,m) + 0.75*this->Distribution(i,j,k,l,m);
                      }

              RKState = 1;
              break;
              
      case 1: forwardDist(this->Distribution, T1Dist, 2*timestep/3.);         
              this->Distribution.swap(T1Dist);

              RKState = 2;
              //VlasovDiagnostic = true;
              break;
              
      case 2: forwardDist(this->Distribution, T2Dist, 3*timestep/4.); 
              this->Distribution.swap(T2Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) 
                      {
                        T1Dist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
                      } 
              RKState = 0;
              break;
  }
}

//===================================================================
//==============FiniteDiff4SchemeRK4 =================================
//===================================================================


extern int myLocalFiniteDiff4SchemeRK3Counter;

template<class ForceField>
void FiniteDiff4SchemeRK3<ForceField>
        ::forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep)
{
    myLocalFiniteDiff4SchemeRK3Counter++;
//    std::cerr << "forwardDist "<<myLocalFiniteDiff4SchemeRK3Counter<<"\n";
    const int *UBound = this->Distribution.getHigh();
    const int *LBound = this->Distribution.getLow();
    
    PositionI Xi;
    VelocityI Vi;
    
    int lx = LBound[0]+2;
    int bx = UBound[0]-2;
    int ly = LBound[1]+2;
    int by = UBound[1]-2;
    int lvx = LBound[2]+1;
    int bvx = UBound[2]-1;
    int lvy = LBound[3]+1;
    int bvy = UBound[3]-1;
    int lvz = LBound[4]+1;
    int bvz = UBound[4]-1;
    
    double dvx = this->deltaVx();
    double dvy = this->deltaVy();
    double dvz = this->deltaVz();
    
//    std::cerr << "forwardDist -- in loop\n";
    for (Vi[2] = lvz; Vi[2] <= bvz; ++Vi[2])
      for (Vi[1] = lvy; Vi[1] <= bvy; ++Vi[1])	
        for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0])
          for (Xi[1] = ly;  Xi[1] <= by;  ++Xi[1])
            for (Xi[0] = lx;  Xi[0] <= bx;  ++Xi[0])
            {
//              if (myLocalFiniteDiff4SchemeRK3Counter==8) {
//                std::cerr << "forwardDist (" << Xi[0]<<","<< Xi[1]<<","<< Vi[0]<<","<< Vi[1]<<","<< Vi[2]<<"): lhs = "
//                  <<lhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) <<"; rhs="
//                  <<rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) << "\n" << flush;
//              }
 
              VelocityD Vel = this->velocity(Vi);
              VelocityD F = this->Force(Xi,Vel,timestep);

              double dfx = (  rhs(Xi[0]-2, Xi[1], Vi[0], Vi[1], Vi[2])
                          - 8*rhs(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])
                          + 8*rhs(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          -   rhs(Xi[0]+2, Xi[1], Vi[0], Vi[1], Vi[2])) / (12.*this->dx[0]);
 
              double dfy = (  rhs(Xi[0], Xi[1]-2, Vi[0], Vi[1], Vi[2])
                          - 8*rhs(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])
                          + 8*rhs(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          -   rhs(Xi[0], Xi[1]+2, Vi[0], Vi[1], Vi[2])) / (12.*this->dx[1]);

              double dfvx;
              double dfvy = 0.01;
              double dfvz = 0.01;
              
              if ((Vi[0]==lvx) || (Vi[0]==bvx))
              {
                dfvx = (rhs(Xi[0], Xi[1], Vi[0]+1, Vi[1], Vi[2])
                      - rhs(Xi[0], Xi[1], Vi[0]-1, Vi[1], Vi[2])) / (2*dvx);
              } else
              {
                dfvx = ( rhs(Xi[0], Xi[1], Vi[0]-2, Vi[1], Vi[2])
                     - 8*rhs(Xi[0], Xi[1], Vi[0]-1, Vi[1], Vi[2])
                     + 8*rhs(Xi[0], Xi[1], Vi[0]+1, Vi[1], Vi[2])
                     -   rhs(Xi[0], Xi[1], Vi[0]+2, Vi[1], Vi[2])) / (12.*dvx);
              }

              if ((Vi[1]==lvy) || (Vi[1]==bvy))
              {
                dfvy = (rhs(Xi[0], Xi[1], Vi[0], Vi[1]+1, Vi[2])
                      - rhs(Xi[0], Xi[1], Vi[0], Vi[1]-1, Vi[2])) / (2*dvy);
              } else
              {
                dfvy = ( rhs(Xi[0], Xi[1], Vi[0], Vi[1]-2, Vi[2])
                     - 8*rhs(Xi[0], Xi[1], Vi[0], Vi[1]-1, Vi[2])
                     + 8*rhs(Xi[0], Xi[1], Vi[0], Vi[1]+1, Vi[2])
                     -   rhs(Xi[0], Xi[1], Vi[0], Vi[1]+2, Vi[2])) / (12.*dvy);

              }

              if ((Vi[2]==lvz) || (Vi[2]==bvz))
              {
                dfvz = (rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+1)
                      - rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-1)) / (2*dvz);
              } else
              {
                dfvz = ( rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-2)
                     - 8*rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-1)
                     + 8*rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+1)
                     -   rhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+2)) / (12.*dvz);
              }
                           
              lhs(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) 
                += - ( timestep*(Vel[0]*dfx + Vel[1]*dfy)
                          + (F[0]*dfvx + F[1]*dfvy + F[2]*dfvz)
                        );
            }
//    std::cerr << "forwardDist -- Boundary\n";
    this->boundary->exchangeX(lhs);
    this->boundary->exchangeY(lhs);

//    std::cerr << "forwardDist -- denErr\n";
    double denErr = this->densityError(lhs)/numDensity;
    this->correctDensityError(denErr, lhs);

//    std::cerr << "forwardDist -- done\n";
}

template<class ForceField>
void FiniteDiff4SchemeRK3<ForceField>::initDist()
{
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
    
  T1Dist.resize(LBound,UBound);
  T2Dist.resize(LBound,UBound);

  for (int i=LBound[0]; i<=UBound[0]; ++i)
    for (int j=LBound[1]; j<=UBound[1]; ++j)
      for (int k=LBound[2]; k<=UBound[2]; ++k)
        for (int l=LBound[3]; l<=UBound[3]; ++l) 
          for (int m=LBound[4]; m<=UBound[4]; ++m)
          {
            T1Dist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
          }
  if (this->boundary->master()) std::cout << "Initializing Density Error\n";
  numDensity = this->densityError(this->Distribution);
}


template<class ForceField>
void FiniteDiff4SchemeRK3<ForceField>
        ::advanceAll(double timestep) 
{
//  std::cerr << "advanceAll "<< RKState <<"\n";
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();

  switch (RKState) {
      case 0: forwardDist(this->Distribution, T1Dist, timestep/3.);         
              this->Distribution.swap(T1Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                      {
                        T2Dist(i,j,k,l,m) = 0.25*T1Dist(i,j,k,l,m) + 0.75*this->Distribution(i,j,k,l,m);
                      }

              RKState = 1;
              break;
              
      case 1: forwardDist(this->Distribution, T1Dist, 2*timestep/3.);         
              this->Distribution.swap(T1Dist);

              RKState = 2;
              //VlasovDiagnostic = true;
              break;
              
      case 2: forwardDist(this->Distribution, T2Dist, 3*timestep/4.); 
              this->Distribution.swap(T2Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) 
                      {
                        T1Dist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
                      } 
              RKState = 0;
              break;
  }
//  std::cerr << "advanceAll -- done\n";

}
