

template<class ForceField>
void FiniteDiffScheme<ForceField>
        ::advanceAll(double timestep) {
    
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
              
              double &D = Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]);
              double dfx = (Distribution(Xi[0]+1, Xi[1], Vi[0], Vi[1], Vi[2])
                          - Distribution(Xi[0]-1, Xi[1], Vi[0], Vi[1], Vi[2])) / (2*dx[0]);

              double dfy = (Distribution(Xi[0], Xi[1]+1, Vi[0], Vi[1], Vi[2])
                          - Distribution(Xi[0], Xi[1]-1, Vi[0], Vi[1], Vi[2])) / (2*dx[1]);

              double dfvx = (Distribution(Xi[0], Xi[1], Vi[0]+1, Vi[1], Vi[2])
                           - Distribution(Xi[0], Xi[1], Vi[0]-1, Vi[1], Vi[2])) / (2*dvx);

              double dfvy = (Distribution(Xi[0], Xi[1], Vi[0], Vi[1]+1, Vi[2])
                           - Distribution(Xi[0], Xi[1], Vi[0], Vi[1]-1, Vi[2])) / (2*dvy);

              double dfvz = (Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]+1)
                           - Distribution(Xi[0], Xi[1], Vi[0], Vi[1], Vi[2]-1)) / (2*dvz);
                           
              D = D - 2*( timestep*(Vel[0]*dfx + Vel[1]*dfy)
                          + (F[0]*dfvx + F[1]*dfvy + F[2]*dfvz)
                        );
            }
    boundary->exchangeX(Distribution);
    boundary->exchangeY(Distribution);
    
}
