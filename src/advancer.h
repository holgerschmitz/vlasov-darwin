// -*- C++ -*-
// $Id$

#include "vlasovbase.h"
#include "scheme.h"
#include "pparabolic.h"

#ifndef ADVANCER_H
#define ADVANCER_H

template<
  class ForceField, 
  template<class> class Scheme = PosFluxCons3rdOrder
>
class SimpleAdvance  : public Scheme<ForceField> {
  private:
      int InitState;
  protected:
 
      /// Advance the distribution function one timestep
      void advance(double timestep);
    
  public:
      SimpleAdvance(SpeciesData &data) 
          : InitState(-1), Scheme<ForceField>(data) 
      {
        if (init->restart()) InitState = 0;
      }
};

template<
  class ForceField, 
  template<class> class Scheme = PosFluxCons3rdOrder
>
class LeapFrogAdvance  : public Scheme<ForceField> {
  private:
      int InitState;
  protected:

      /** @brief Half of the distribution advance.
       *  Only advances the spatial direction for half a timestep
       */
      void advanceStepA(double timestep);

      /** @brief Half of the distribution advance.
       *  Only advances the velocity directions one timestep
       *  and then the spatial direction for half a timestep.
       */
      void advanceStepB(double timestep);

      /** @brief The full distribution advance.
       *  Only advances the velocity directions one timestep
       *  and then the spatial direction for a timestep.
       */
      void advanceStepFull(double timestep);
      /// Advance the distribution function one timestep
      void advance(double timestep);
    
  public:
      LeapFrogAdvance(SpeciesData &data) 
          : InitState(-2), Scheme<ForceField>(data) 
      {
        if (init->restart()) InitState = 0;
      }
};

template<
  class ForceField, 
  template<class> class Scheme = PosFluxCons3rdOrder
>
class SimpleLeapFrogAdvanceBase  : public Scheme<ForceField> {
  protected:

      /** @brief Half of the distribution advance.
       *  Only advances the spatial direction for half a timestep
       */
      void advanceStepA(double timestep);

      /** @brief Half of the distribution advance.
       *  Only advances the velocity directions one timestep
       *  and then the spatial direction for half a timestep.
       */
      void advanceStepB(double timestep);

      /** @brief The full distribution advance.
       *  Only advances the velocity directions one timestep
       *  and then the spatial direction for a timestep.
       */
      void advanceStepFull(double timestep);
    
  public:
      SimpleLeapFrogAdvanceBase(SpeciesData &data) : Scheme<ForceField>(data) {}
};

template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class SimpleLeapFrogAdvance : public SimpleLeapFrogAdvanceBase<ForceField,Scheme> {
  private:
      int InitState;
  public:
    SimpleLeapFrogAdvance(SpeciesData &data) 
      : InitState(-2), 
        SimpleLeapFrogAdvanceBase<ForceField,Scheme>(data) 
    {
      if (init->restart()) InitState = 0;
    }
    
    /// Advance the distribution function one timestep
    void advance(double timestep);
};

template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class RungeKuttaAdvance : public SimpleLeapFrogAdvanceBase<ForceField,Scheme> {
  private:
      /** The state of the Runge Kutta integrator.
       *  Also used by the leap-frog integrator for the start-up phase
       */
      int RKState;
      /// The Runge Kutta Temporaries
      VlasovDist T1, T2;
      /// Another temporary needed for the Runge Kutta integration
      VlasovDist TempDist;

  public:
      RungeKuttaAdvance(SpeciesData &data) 
          : RKState(-2), SimpleLeapFrogAdvanceBase<ForceField,Scheme>(data) {}
      void initializeAdvancer();
      /// Advance the distribution function one timestep
      void advance(double timestep);
};

template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class RungeKuttaBAdvance : public SimpleLeapFrogAdvanceBase<ForceField,Scheme> {
  private:
      /** The state of the Runge Kutta integrator.
       *  Also used by the leap-frog integrator for the start-up phase
       */
      int RKState;
      /// The Runge Kutta Temporaries
      VlasovDist T1, T2;
      /// Another temporary needed for the Runge Kutta integration
      VlasovDist TempDist;

  public:
      RungeKuttaBAdvance(SpeciesData &data) 
          : RKState(-2), SimpleLeapFrogAdvanceBase<ForceField,Scheme>(data) {}
      void initializeAdvancer();
      /// Advance the distribution function one timestep
      void advance(double timestep);
};

#include "advancer_temp.cpp"

#endif
