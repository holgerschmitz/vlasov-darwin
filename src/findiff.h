// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
#include "vlasovbase.h"
//-----------------------------------------------------------------------------
#ifndef FINDIFF_H
#define FINDIFF_H
//-----------------------------------------------------------------------------
  /** @file findiff.h
	* @brief implements various finite difference schemes
	*
	* since finite difference schemes to not work with the advancer class
      * explicit Runge Kutta integrators are used for time advancement
	*/
///@todo implement an independent RK integrator for plugging into the finite difference schemes
//-----------------------------------------------------------------------------
//FiniteDiffSchemeRK
///class for implementing the finite difference scheme, derived from forcefield
template<class ForceField>
class FiniteDiffScheme : public ForceField {
  protected:
      /// A maximum value of the distribution function
      double f_infty;
      ///temporary distribution function
      VlasovDist TempDist;
	///current distribution function
      VlasovDist CurrentDist;
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      ///initializer
      void initDist();
  public:
      ///constructor
      FiniteDiffScheme(SpeciesData &data) : ForceField(data) {}
};
//FiniteDiffScheme
//-----------------------------------------------------------------------------
//FiniteDiffSchemeRK2
/** @brief Finite difference scheme,
  *
  * using a Runge-Kutta integrator of 2nd order
  */
template<class ForceField>
class FiniteDiffSchemeRK2 : public ForceField {
  private:
	///advance the distribution, to be used by advanceAll
      void forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep);
	///state of the integrator
      int RKState;
  protected:
      /// A maximum value of the distribution function
      double f_infty;
	///temporary distribution function
      VlasovDist D1;
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      ///initializer
      void initDist();
  public:
  ///constructor
      FiniteDiffSchemeRK2(SpeciesData &data)
        : ForceField(data), RKState(0) {}
};
//FiniteDiffSchemeRK2
//-----------------------------------------------------------------------------
//FiniteDiffSchemeRK3
/** @brief Finite difference scheme,
 *
 * using a Runge-Kutta integrator of 3rd order
 */
template<class ForceField>
class FiniteDiffSchemeRK3 : public ForceField {
  private:
	///advance one direction
      void forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep);
	///state of the integrator
      int RKState;
	///numerical density
      double numDensity;
  protected:
      /// A maximum value of the distribution function
      double f_infty;
	VlasovDist T1Dist, ///< temporary distribution function
		     T2Dist; ///< temporary distribution function
      
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      ///initializer
      void initDist();
  public:
  ///constructor
      FiniteDiffSchemeRK3(SpeciesData &data)
        : ForceField(data), RKState(0) {}
};
//FiniteDiffSchemeRK3
//-----------------------------------------------------------------------------
//FiniteDiff4SchemeRK3
/** @brief finite difference scheme,
 *
 * using a Runge-Kutta integrator of 4th order
 */
template<class ForceField>
class FiniteDiff4SchemeRK3 : public ForceField {
  private:
	///advnance one direction
      void forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep);
	///stae of the integrator
      int RKState;
	///numerical Density
      double numDensity;
  protected:
      /// A maximum value of the distribution function
      double f_infty;
      
	VlasovDist T1Dist, ///<temporary distribution function
	T2Dist; ///<temporary distribution function
      
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      ///initializer
      void initDist();
  public:
  ///constructor
      FiniteDiff4SchemeRK3(SpeciesData &data)
        : ForceField(data), RKState(0) {}
};
//FiniteDiff4SchemeRK3
//-----------------------------------------------------------------------------
//include template implementations
#include "findiff.t"
//-----------------------------------------------------------------------------
#endif //FINDIFF_H
