template<
  class ForceField, 
  template<class> class Scheme
>
void LeapFrogAdvanceBase<ForceField,Scheme>
        ::advanceStepA(double timestep) {
        
//    InterpolationInitStep(Distribution);
    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
}

template<
  class ForceField, 
  template<class> class Scheme
>
void LeapFrogAdvanceBase<ForceField,Scheme>
        ::advanceStepB(double timestep) {

    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_z(timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);

    advanceSpace_x(0.5*timestep);
    advanceSpace_y(0.5*timestep);
    
}

template<
  class ForceField, 
  template<class> class Scheme
>
void LeapFrogAdvanceBase<ForceField,Scheme>
        ::advanceStepFull(double timestep) {

    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_z(timestep);
    advanceVel_x(0.25*timestep);
    advanceVel_y(0.5*timestep);
    advanceVel_x(0.25*timestep);

    advanceSpace_x(timestep);
    advanceSpace_y(timestep);
    
}

template<
  class ForceField, 
  template<class> class Scheme
>
void LeapFrogAdvance<ForceField,Scheme>::advance(double timestep) {
  switch (InitState) {
      case -2: InitState = -1;
               break;
      case -1: InitState = 0;
              advanceStepA(timestep); 
              break;
      case 0: advanceStepFull(timestep); 
              break;
  }
}


template<
  class ForceField, 
  template<class> class Scheme
>
void RungeKuttaAdvance<ForceField,Scheme>::initializeAdvancer() {
  PhasePositionI low = boundary->DistLow();
  PhasePositionI high = boundary->DistHigh();
  
  TempDist.resize(low.Data(),high.Data());
  T1.resize(low.Data(),high.Data());
  T2.resize(low.Data(),high.Data());
}


template<
  class ForceField, 
  template<class> class Scheme
>
void RungeKuttaAdvance<ForceField,Scheme>::advance(double timestep) {
  const int *UBound = Distribution.getHigh();
  const int *LBound = Distribution.getLow();

  switch (RKState) {
      case -2: RKState = -1;
              initializeAdvancer();
              break;
      case -1: RKState = 0;
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                        TempDist(i,j,k,l,m) = Distribution(i,j,k,l,m);
              advanceStepA(timestep/3.); 
              break;
      case 0: advanceStepB(timestep/3.);         // Distribution = c0 = d1

              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                        T1(i,j,k,l,m) = Distribution(i,j,k,l,m);


              advanceStepA(2*timestep/3.); 
              RKState = 1;
              break;
      case 1: advanceStepB(2*timestep/3.);        // Distribution = c1

              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) {
                        Distribution(i,j,k,l,m) 
                          = Distribution(i,j,k,l,m) + TempDist(i,j,k,l,m) - T1(i,j,k,l,m);
                        T2(i,j,k,l,m) = Distribution(i,j,k,l,m); // = d2
                      }
              advanceStepA(3*timestep/4.);

              RKState = 2;
              break;
      case 2: advanceStepB(3*timestep/4.);      // Distribution = c2
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) {
                        Distribution(i,j,k,l,m) 
                            =  (1/4.)*TempDist(i,j,k,l,m) 
                             + (3/4.)*T1(i,j,k,l,m) 
                             +  Distribution(i,j,k,l,m)
                             - T2(i,j,k,l,m); 
                        TempDist(i,j,k,l,m) = Distribution(i,j,k,l,m);
                      } 
              advanceStepA(timestep/3.); 
              RKState = 0;
              break;
  }
}

