// -*- C++ -*-
// $Id$

#include "vlasovbase.h"
#include "scheme.h"

#ifndef ADVANCER_H
#define ADVANCER_H


template<
  class ForceField, 
  template<class> class Scheme = PosFluxCons3rdOrder
>
class SimpleLeapFrogAdvance  : public Scheme<ForceField> {
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
      SimpleLeapFrogAdvance(SpeciesData &data) 
          : InitState(-2), Scheme<ForceField>(data) {}
};

template<
  class ForceField, 
  template<class> class Scheme = PosFluxCons3rdOrder
>
class LeapFrogAdvanceBase  : public Scheme<ForceField> {
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
      LeapFrogAdvanceBase(SpeciesData &data) : Scheme<ForceField>(data) {}
};

template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class LeapFrogAdvance : public LeapFrogAdvanceBase<ForceField,Scheme> {
  private:
      int InitState;
  public:
      LeapFrogAdvance(SpeciesData &data) 
          : InitState(-2), LeapFrogAdvanceBase<ForceField,Scheme>(data) {}
      /// Advance the distribution function one timestep
      void advance(double timestep);
};

template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class RungeKuttaAdvance : public LeapFrogAdvanceBase<ForceField,Scheme> {
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
          : RKState(-2), LeapFrogAdvanceBase<ForceField,Scheme>(data) {}
      void initializeAdvancer();
      /// Advance the distribution function one timestep
      void advance(double timestep);
};

#include "advancer_temp.cpp"

#endif
