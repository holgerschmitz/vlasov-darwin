
template<
  class ForceField, 
  template<class> class Scheme
>
void SimpleAdvance<ForceField,Scheme>::advance(double timestep) {
  switch (InitState) {
      case -1: InitState = 0;
               this->initDist();
               break;
      case  0: this->advanceAll(timestep);
               break;
  }
}


template<
  class ForceField, 
  template<class> class Scheme
>
void LeapFrogAdvance<ForceField,Scheme>
        ::advanceStepA(double timestep) {
        
    this->advanceSpace_x(0.5*timestep);
    this->advanceSpace_y(0.5*timestep);
    
}

template<
  class ForceField, 
  template<class> class Scheme
>
void LeapFrogAdvance<ForceField,Scheme>
        ::advanceStepB(double timestep) {

    this->advanceVel_x(0.25*timestep);
    this->advanceVel_y(0.5*timestep);
    this->advanceVel_x(0.25*timestep);
    this->advanceVel_z(timestep);
    this->advanceVel_x(0.25*timestep);
    this->advanceVel_y(0.5*timestep);
    this->advanceVel_x(0.25*timestep);

    this->advanceSpace_x(0.5*timestep);
    this->advanceSpace_y(0.5*timestep);
    
}

template<
  class ForceField, 
  template<class> class Scheme
>
void LeapFrogAdvance<ForceField,Scheme>
        ::advanceStepFull(double timestep) {

    this->advanceVel_x(0.25*timestep);
    this->advanceVel_y(0.5*timestep);
    this->advanceVel_x(0.25*timestep);
    this->advanceVel_z(timestep);
    this->advanceVel_x(0.25*timestep);
    this->advanceVel_y(0.5*timestep);
    this->advanceVel_x(0.25*timestep);

    this->advanceSpace_x(timestep);
    this->advanceSpace_y(timestep);
    
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
              this->advanceStepA(timestep); 
              break;
      case 0: this->advanceStepFull(timestep); 
              break;
  }
  
}


template<
  class ForceField, 
  template<class> class Scheme
>
void SimpleLeapFrogAdvanceBase<ForceField,Scheme>
        ::advanceStepA(double timestep) {
//    std::cerr << "advanceStepA x\n";
    this->advanceSpace_x(0.5*timestep);
//    std::cerr << "advanceStepA y\n";
    this->advanceSpace_y(0.5*timestep);
//    std::cerr << "advanceStepA done\n";
}

template<
  class ForceField, 
  template<class> class Scheme
>
void SimpleLeapFrogAdvanceBase<ForceField,Scheme>
        ::advanceStepB(double timestep) {

    this->advanceVel_z(timestep);
    this->advanceVel_y(timestep);
    this->advanceVel_x(timestep);

    this->advanceSpace_x(0.5*timestep);
    this->advanceSpace_y(0.5*timestep);
}

template<
  class ForceField, 
  template<class> class Scheme
>
void SimpleLeapFrogAdvanceBase<ForceField,Scheme>
        ::advanceStepFull(double timestep) {

    this->advanceVel_z(timestep);
    this->advanceVel_y(timestep);
    this->advanceVel_x(timestep);

//    std::cerr << "advanceStepFull x\n";
    this->advanceSpace_x(timestep);
//    std::cerr << "advanceStepFull y\n";
    this->advanceSpace_y(timestep);
//    std::cerr << "advanceStepFull done\n";
}

template<
  class ForceField, 
  template<class> class Scheme
>
void SimpleLeapFrogAdvance<ForceField,Scheme>::advance(double timestep) {
  switch (InitState) {
      case -2: InitState = -1;
               break;
               
      case -1: InitState = 0;
               if (restarting) this->advanceStepFull(timestep);
               else this->advanceStepA(timestep); 
               Parameters::instance().setRestart(false);
               break;
               
      case 0:  this->advanceStepFull(timestep);
               break;
  }
}


template<
  class ForceField, 
  template<class> class Scheme
>
void RungeKuttaAdvance<ForceField,Scheme>::initializeAdvancer() {
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
    
  TempDist.resize(LBound,UBound);
  T1.resize(LBound,UBound);
  T2.resize(LBound,UBound);
}


template<
  class ForceField, 
  template<class> class Scheme
>
void RungeKuttaAdvance<ForceField,Scheme>::advance(double timestep) {
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();

  switch (RKState) {
      case -2: RKState = -1;
              this->initializeAdvancer();
              break;
      case -1: RKState = 0;
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                        TempDist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
              this->advanceStepA(timestep/3.); 
              break;
      case 0: this->advanceStepB(timestep/3.);         // Distribution = c0 = d1

              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                        T1(i,j,k,l,m) = this->Distribution(i,j,k,l,m);


              this->advanceStepA(2*timestep/3.); 
              RKState = 1;
              break;
      case 1: this->advanceStepB(2*timestep/3.);        // Distribution = c1

              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) {
                        this->Distribution(i,j,k,l,m) 
                          = this->Distribution(i,j,k,l,m) + TempDist(i,j,k,l,m) - T1(i,j,k,l,m);
                        T2(i,j,k,l,m) = this->Distribution(i,j,k,l,m); // = d2
                      }
              this->advanceStepA(3*timestep/4.);

              RKState = 2;
              //VlasovDiagnostic = true;
              break;
      case 2: this->advanceStepB(3*timestep/4.);      // Distribution = c2
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) {
                        this->Distribution(i,j,k,l,m) 
                            =  (1/4.)*TempDist(i,j,k,l,m) 
                             + (3/4.)*T1(i,j,k,l,m) 
                             +  this->Distribution(i,j,k,l,m)
                             - T2(i,j,k,l,m); 
                        TempDist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
                      } 
              this->advanceStepA(timestep/3.); 
              RKState = 0;
              break;
  }
}

template<
  class ForceField, 
  template<class> class Scheme
>
void RungeKuttaBAdvance<ForceField,Scheme>::initializeAdvancer() {
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();
    
  TempDist.resize(LBound,UBound);
  T1.resize(LBound,UBound);
  T2.resize(LBound,UBound);
}


template<
  class ForceField, 
  template<class> class Scheme
>
void RungeKuttaBAdvance<ForceField,Scheme>::advance(double timestep) {
  const int *UBound = this->Distribution.getHigh();
  const int *LBound = this->Distribution.getLow();

  switch (RKState) {
      case -2: RKState = -1;
              this->initializeAdvancer();
              break;
      case -1: RKState = 0;
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                        TempDist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
              this->advanceStepA(timestep/2.); 
              break;
      case 0: this->advanceStepB(timestep/2.);         // Distribution = c0 = d1

              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m)
                        T1(i,j,k,l,m) = this->Distribution(i,j,k,l,m);


              this->advanceStepA(2*timestep/3.); 
              RKState = 1;
              break;
      case 1: this->advanceStepB(2*timestep/3.);        // Distribution = c1

              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) {
                        this->Distribution(i,j,k,l,m) 
                          = fabs(this->Distribution(i,j,k,l,m) + TempDist(i,j,k,l,m) - T1(i,j,k,l,m));
                        T2(i,j,k,l,m) = this->Distribution(i,j,k,l,m); // = d2
                      }
              this->advanceStepA(3*timestep/4.);

              RKState = 2;
              //VlasovDiagnostic = true;
              break;
      case 2: this->advanceStepB(3*timestep/4.);      // Distribution = c2
              for (int i=LBound[0]; i<=UBound[0]; ++i)
                for (int j=LBound[1]; j<=UBound[1]; ++j)
                  for (int k=LBound[2]; k<=UBound[2]; ++k)
                    for (int l=LBound[3]; l<=UBound[3]; ++l) 
                      for (int m=LBound[4]; m<=UBound[4]; ++m) {
                        this->Distribution(i,j,k,l,m) 
                            =  fabs((1/4.)*TempDist(i,j,k,l,m) 
                             + (3/4.)*T1(i,j,k,l,m) 
                             +  this->Distribution(i,j,k,l,m)
                             - T2(i,j,k,l,m)); 
                        TempDist(i,j,k,l,m) = this->Distribution(i,j,k,l,m);
                      } 
              this->advanceStepA(timestep/3.); 
              RKState = 0;
              break;
  }
}

