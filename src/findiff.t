//================================================================
//========       FiniteDiffScheme        =========================
//================================================================

template<class ForceField>
void FiniteDiffScheme<ForceField>::initDist()
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();
    
  TempDist.resize(LBound,UBound);
  CurrentDist.resize(LBound,UBound);
}

template<class ForceField>
void FiniteDiffScheme<ForceField>
        ::advanceAll(double timestep) {
    
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();

                                       //  T | C | D |
    Distribution.swap(CurrentDist);    //  T | D | C |
    Distribution.swap(TempDist);       //  D | T | C |
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
    
    double dvx = deltaVx();
    double dvy = deltaVy();
    double dvz = deltaVz();
    
    for (Vi[2] = lvz; Vi[2] <= bvz; ++Vi[2])
      for (Vi[1] = lvy; Vi[1] <= bvy; ++Vi[1])	
        for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0])
          for (Xi[1] = ly;  Xi[1] <= by;  ++Xi[1])
            for (Xi[0] = lx;  Xi[0] <= bx;  ++Xi[0])
            {
              VelocityD Vel = velocity(Vi);
              VelocityD F = Force(Xi,Vel,timestep);
              
              double dfx = (CurrentDist(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          - CurrentDist(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])) / (2*dx[0]);

              double dfy = (CurrentDist(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          - CurrentDist(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])) / (2*dx[1]);

              double dfvx = (CurrentDist(Xi[0], Xi[1], Vi[0]+1, Vi[1], Vi[2])
                           - CurrentDist(Xi[0], Xi[1], Vi[0]-1, Vi[1], Vi[2])) / (2*dvx);

              double dfvy = (CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1]+1, Vi[2])
                           - CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1]-1, Vi[2])) / (2*dvy);

              double dfvz = (CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+1)
                           - CurrentDist(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-1)) / (2*dvz);
                           
              Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]) 
                = TempDist(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2])
                   - 2*( timestep*(Vel[0]*dfx + Vel[1]*dfy)
                          + (F[0]*dfvx + F[1]*dfvy + F[2]*dfvz)
                        );
            }
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);
    
}

//===================================================================
//==============FiniteDiffSchemeRK2 =================================
//===================================================================

template<class ForceField>
void FiniteDiffSchemeRK2<ForceField>::initDist()
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();
    
  D1.resize(LBound,UBound);
}

template<class ForceField>
void FiniteDiffSchemeRK2<ForceField>
        ::forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep)
{
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    
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
    
    double dvx = deltaVx();
    double dvy = deltaVy();
    double dvz = deltaVz();
    
    for (Vi[2] = lvz; Vi[2] <= bvz; ++Vi[2])
      for (Vi[1] = lvy; Vi[1] <= bvy; ++Vi[1])	
        for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0])
          for (Xi[1] = ly;  Xi[1] <= by;  ++Xi[1])
            for (Xi[0] = lx;  Xi[0] <= bx;  ++Xi[0])
            {
              VelocityD Vel = velocity(Vi);
              VelocityD F = Force(Xi,Vel,timestep);
              
              double dfx = (rhs(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])) / (2*dx[0]);

              double dfy = (rhs(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])) / (2*dx[1]);

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
    boundary->exchangeX(lhs);
    boundary->exchangeY(lhs);
}

template<class ForceField>
void FiniteDiffSchemeRK2<ForceField>
        ::advanceAll(double timestep) 
{   
  switch (RKState)
  {
    case 0:
      forwardDist(D1, Distribution, 0.5*timestep);
      RKState = 1;
      break;
    case 1:
      forwardDist(Distribution, D1, timestep);
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
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    
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
    
    double dvx = deltaVx();
    double dvy = deltaVy();
    double dvz = deltaVz();
    
    for (Vi[2] = lvz; Vi[2] <= bvz; ++Vi[2])
      for (Vi[1] = lvy; Vi[1] <= bvy; ++Vi[1])	
        for (Vi[0] = lvx; Vi[0] <= bvx; ++Vi[0])
          for (Xi[1] = ly;  Xi[1] <= by;  ++Xi[1])
            for (Xi[0] = lx;  Xi[0] <= bx;  ++Xi[0])
            {
              VelocityD Vel = velocity(Vi);
              VelocityD F = Force(Xi,Vel,timestep);
              
              double dfx = (rhs(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])) / (2*dx[0]);

              double dfy = (rhs(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          - rhs(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])) / (2*dx[1]);

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
    boundary->exchangeX(lhs);
    boundary->exchangeY(lhs);
    
    double denErr = densityError(lhs)/numDensity;
    correctDensityError(denErr, lhs);
}

template<class ForceField>
void FiniteDiffSchemeRK3<ForceField>::initDist()
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();
    
  T1Dist.resize(LBound,UBound);
  T2Dist.resize(LBound,UBound);

  for (int i=LBound[0]; i<=UBound[0]; ++i)
    for (int j=LBound[1]; j<=UBound[1]; ++j)
      for (int k=LBound[2]; k<=UBound[2]; ++k)
        for (int l=LBound[3]; l<=UBound[3]; ++l) 
          for (int m=LBound[4]; m<=UBound[4]; ++m)
          {
            T1Dist(i,j,k,l,m) = Distribution(i,j,k,l,m);
          }
  if (boundary->master()) std::cout << "Initializing Density Error\n";
  numDensity = densityError(Distribution);
}


template<class ForceField>
void FiniteDiffSchemeRK3<ForceField>
        ::advanceAll(double timestep) 
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  switch (RKState) {
      case 0: forwardDist(Distribution, T1Dist, timestep/3.);         
              Distribution.swap(T1Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                      {
                        T2Dist(i,j,k,l,m) = 0.25*T1Dist(i,j,k,l,m) + 0.75*Distribution(i,j,k,l,m);
                      }

              RKState = 1;
              break;
              
      case 1: forwardDist(Distribution, T1Dist, 2*timestep/3.);         
              Distribution.swap(T1Dist);

              RKState = 2;
              VlasovDiagnostic = true;
              break;
              
      case 2: forwardDist(Distribution, T2Dist, 3*timestep/4.); 
              Distribution.swap(T2Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) 
                      {
                        T1Dist(i,j,k,l,m) = Distribution(i,j,k,l,m);
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
    const int *UBound = Distribution.getHigh();
    const int *LBound = Distribution.getLow();
    
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
    
    double dvx = deltaVx();
    double dvy = deltaVy();
    double dvz = deltaVz();
    
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
 
              VelocityD Vel = velocity(Vi);
              VelocityD F = Force(Xi,Vel,timestep);

              double dfx = (  rhs(Xi[0]-2, Xi[1], Vi[0], Vi[1], Vi[2])
                          - 8*rhs(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])
                          + 8*rhs(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          -   rhs(Xi[0]+2, Xi[1], Vi[0], Vi[1], Vi[2])) / (12.*dx[0]);
 
              double dfy = (  rhs(Xi[0], Xi[1]-2, Vi[0], Vi[1], Vi[2])
                          - 8*rhs(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])
                          + 8*rhs(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          -   rhs(Xi[0], Xi[1]+2, Vi[0], Vi[1], Vi[2])) / (12.*dx[1]);

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
    boundary->exchangeX(lhs);
    boundary->exchangeY(lhs);

//    std::cerr << "forwardDist -- denErr\n";
    double denErr = densityError(lhs)/numDensity;
    correctDensityError(denErr, lhs);

//    std::cerr << "forwardDist -- done\n";
}

template<class ForceField>
void FiniteDiff4SchemeRK3<ForceField>::initDist()
{
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();
    
  T1Dist.resize(LBound,UBound);
  T2Dist.resize(LBound,UBound);

  for (int i=LBound[0]; i<=UBound[0]; ++i)
    for (int j=LBound[1]; j<=UBound[1]; ++j)
      for (int k=LBound[2]; k<=UBound[2]; ++k)
        for (int l=LBound[3]; l<=UBound[3]; ++l) 
          for (int m=LBound[4]; m<=UBound[4]; ++m)
          {
            T1Dist(i,j,k,l,m) = Distribution(i,j,k,l,m);
          }
  if (boundary->master()) std::cout << "Initializing Density Error\n";
  numDensity = densityError(Distribution);
}


template<class ForceField>
void FiniteDiff4SchemeRK3<ForceField>
        ::advanceAll(double timestep) 
{
//  std::cerr << "advanceAll "<< RKState <<"\n";
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  switch (RKState) {
      case 0: forwardDist(Distribution, T1Dist, timestep/3.);         
              Distribution.swap(T1Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                      {
                        T2Dist(i,j,k,l,m) = 0.25*T1Dist(i,j,k,l,m) + 0.75*Distribution(i,j,k,l,m);
                      }

              RKState = 1;
              break;
              
      case 1: forwardDist(Distribution, T1Dist, 2*timestep/3.);         
              Distribution.swap(T1Dist);

              RKState = 2;
              VlasovDiagnostic = true;
              break;
              
      case 2: forwardDist(Distribution, T2Dist, 3*timestep/4.); 
              Distribution.swap(T2Dist);
              
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) 
                      {
                        T1Dist(i,j,k,l,m) = Distribution(i,j,k,l,m);
                      } 
              RKState = 0;
              break;
  }
//  std::cerr << "advanceAll -- done\n";

}
