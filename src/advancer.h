// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
#include "vlasovbase.h"
#include "scheme.h"
#include "pparabolic.h"
#include "globals.h"
//-----------------------------------------------------------------------------
/** @file advancer.h
  * @brief implements various advancers
  *
  * The Advancer calls the used scheme for spatial integration, so splitting the 
  * time step to carry out a spatial integration is possible.
  *  Available schemes are Runge Kutta, Euler, leap frog
  */

//-----------------------------------------------------------------------------
#ifndef ADVANCER_H
#define ADVANCER_H
//-----------------------------------------------------------------------------
//SimpleAdvance	
///simple advancer template, Euler
template<
  class ForceField, 
  template<class> class Scheme = PosFluxCons3rdOrder
>
class SimpleAdvance  : public Scheme<ForceField> {
  private:
     /** The state of the  integrator.
	 */
      int InitState;
  protected:
 
      /// Advance the distribution function one timestep
      void advance(double timestep);
    
  public:
      ///constructor
      SimpleAdvance(SpeciesData &data) 
          : InitState(-1), Scheme<ForceField>(data) 
      {
        //if (init->restart()) InitState = 0;
      }
};
//SimpleAdvance
//-----------------------------------------------------------------------------
//LeapFrogAdvance

///leap frog advancer
template<
  class ForceField, 
  template<class> class Scheme = PosFluxCons3rdOrder
>
class LeapFrogAdvance  : public Scheme<ForceField> {
  private:
    /** The state of the integrator.
      * Used for the start-up phase
	*/
      int InitState;
  protected:

     /** @brief Half of the distribution advance.
	 *
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
      ///constructor
      LeapFrogAdvance(SpeciesData &data) 
          : InitState(-2), Scheme<ForceField>(data) 
      {
        if (this->init->restart()) InitState = 0;
      }
};
//LeapFrogAdvance
//-----------------------------------------------------------------------------
//SimpleLeapFrogAdvanceBase

///simple leap frog advancer, base class
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
       *
       *  Only advances the velocity directions one timestep
       *  and then the spatial direction for a timestep.
       */
      void advanceStepFull(double timestep);
    
  public:
  	///constructor, takes species data for initializing the Scheme
      SimpleLeapFrogAdvanceBase(SpeciesData &data) : Scheme<ForceField>(data) {}
};
//SimpleLeapFrogBase
//-----------------------------------------------------------------------------
//SimpleLeapFrog

/** @brief simple leap frog advancer 
  *
  * here the splitting is just in velocity and spatial inegration
  */
template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class SimpleLeapFrogAdvance : public SimpleLeapFrogAdvanceBase<ForceField,Scheme> {
  private:
    /** The state of the Runge Kutta integrator.
      *  Also used by the leap-frog integrator for the start-up phase
	*/
      int InitState;
	/// "true" if its a restart
      bool restarting;
  public:
  ///constructor
    SimpleLeapFrogAdvance(SpeciesData &data) 
      : InitState(-2), 
        SimpleLeapFrogAdvanceBase<ForceField,Scheme>(data) 
    {
      //if (this->init->restart()) InitState = 0;
      restarting = (this->init->restart());
      if (restarting) std::cout << "RESTARTING" << std::endl;
    }
    
    /// Advance the distribution function one timestep
    void advance(double timestep);
};
//SimpleLeapFrog
//-----------------------------------------------------------------------------
//RungeKuttaAdvance

/** @brief First Runge Kutta advancer, derived from SimpleLeapFrog
 */
template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class RungeKuttaAdvance : public SimpleLeapFrogAdvanceBase<ForceField,Scheme> {
  private:
      /** The state of the Runge Kutta integrator.
       *  Also used by the leap-frog integrator for the start-up phase
       */
      int RKState;
      
	VlasovDist T1, ///< First Runge Kutta Temporary
		     T2; ///< Second Runge Kutta Temporary
      /// Another temporary needed for the Runge Kutta integration
      VlasovDist TempDist;

  public:
      ///constructor
      RungeKuttaAdvance(SpeciesData &data) 
          : RKState(-2), SimpleLeapFrogAdvanceBase<ForceField,Scheme>(data) {}
      ///initializer
      void initializeAdvancer();
      /// Advance the distribution function one timestep
      void advance(double timestep);
};
//RungeKuttaAdvance
//-----------------------------------------------------------------------------
//RungeKuttaBAdvance

/** @brief Second Runge Kutta advancer, derived from SimpleLeapFrog
  *
  *
  */
template<class ForceField, template<class> class Scheme = PosFluxCons3rdOrder>
class RungeKuttaBAdvance : public SimpleLeapFrogAdvanceBase<ForceField,Scheme> {
  private:
      /** The state of the Runge Kutta integrator.
       *  Also used by the leap-frog integrator for the start-up phase
       */
      int RKState;
	VlasovDist T1, ///< First Runge Kutta Temporary
		     T2; ///< Second Runge Kutta Temporary
      /// Another temporary needed for the Runge Kutta integration
      VlasovDist TempDist;

  public:
      ///constructor
      RungeKuttaBAdvance(SpeciesData &data) 
          : RKState(-2), SimpleLeapFrogAdvanceBase<ForceField,Scheme>(data) {}
      ///initializer
      void initializeAdvancer();
      /// Advance the distribution function one timestep
      void advance(double timestep);
};
//RungeKuttaBAdvance
//-----------------------------------------------------------------------------
//including the member definitions
#include "advancer_temp.cpp"

#endif //ADVANCER_H
